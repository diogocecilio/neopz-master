# EmbankmentConsolidation — Sect. 6.5, Figs. 10–12 and Table 7

Embankment loading on a Modified Cam-Clay foundation (FLAC3D example *Embankment loading on a Cam-Clay
foundation*): a 50 kPa strip load on a 10 m saturated clay layer, applied in undrained increments and
followed by consolidation up to t = 10⁸ s. This is the function `aterro()` of `gen_data.py`; the elastic
variant is `aterro_elastic.py`. The example also prints the embankment results of Sect. 6.6 (global Newton
iterations).

* Domain: half of the problem, [0, 20] × [0, 10] m, plane strain, 20 × 10 Q8–Q4 elements (serendipity
  displacement, linear pore pressure; 661 displacement nodes, 231 pore pressure nodes, 1553 equations),
  3 × 3 Gauss points (`SetIntegrationOrder(4)`), Fig. 10.
* Material (Table 1): M = 0.888, λ = 0.161, κ = 0.062, v_λ = 2.858, porous elasticity with ν = 0.3,
  uniform p'c0 = 160 kPa. γ_sat = 23 kN/m³, γ_w = 10 kN/m³, K_f = 2·10⁵ kPa, n = 0.3 (α_B = 1,
  1/M_B = n/K_f), mobility k = 10⁻⁹ m²/(kPa s); body force (0, −23) and fluid weight (0, −10).
* Initial state at the depth d = 10 − y: p_w = 10 d, σ'_v = −13 d, σ'_h = 0.7(−23 d) + 10 d = −6.1 d,
  p'c = 160 kPa and v0 = v_λ − λ ln p'c0 + κ ln(p'c0/p'0) at each integration point
  (`TPZMatPoroElastoPlasticUP::InitializeMemory`); hydrostatic nodal pore pressures (`mcc::SetInitialPressure`).
* Boundary conditions: u_x = 0 at x = 0 and x = 20 m (`EDirichletUDirectional`), fixed and impermeable base
  (`EDirichletU`), drained top (`EDirichletP` on coincident line elements over the whole top) and the strip
  load q = 50 kPa on 0 ≤ x ≤ 4 m (`ENeumannU`, scaled by the load factor).
* Loading: 10 undrained increments of the load factor (Δt = 0), then 25 consolidation steps,
  t = 10^(2 + j/4) s, j = 0…24; no predictor, no bisection needed.
* Monitoring: settlements of the top at x = 0, 2, 4, 6 m and the pore pressures pp1 and pp2 (means of the
  four vertices of the elements centred at (0.5, 9.5) and (1.5, 7.5) m).

Three models are solved (`EmbankmentConsolidation::EVariant`):

| Variant | Python | What it is |
|---|---|---|
| `ECamClay` | `aterro()` | the model above: undrained loading and consolidation |
| `ETransposed` | end of `aterro()` (`logUt`) | undrained loading only, with the transposed tangent Dᵀ |
| `EElastic` | `aterro_elastic.py` | p_c = 10⁷ kPa in the state (no yielding), v0 still computed with p'c0 = 160 kPa; gives the plastic share of the final settlement |

The class `EmbankmentConsolidation` follows the structure of the NeoPZ examples:

* `CreateGeoMesh` builds the mesh and its boundary lines with `mcc::CreateRectangleMesh`.
* `CreateCompMesh` builds the displacement, pore pressure and multiphysics meshes. It also sets the
  `TPZMatPoroElastoPlasticUP` material, the boundary conditions and the initial state.
* `Run` holds the analysis, `TPZPoroElastoPlasticUPAnalysis` with `TPZSkylineNSymStructMatrix` and
  `TPZStepSolver` LU. It applies the increments and does the post-processing.
* `RunAll` runs the three models and prints the comparison.

## Running

```
./EmbankmentConsolidation          # CSV files and the VTK series of every converged state
./EmbankmentConsolidation novtk    # without the VTK series
```

Without the VTK series the run takes about 34 s: 15 s for the Cam-Clay model, 5 s for the transposed-tangent
loading and 14 s for the elastic model. The Python script takes 2–3 min. Almost all of the time goes into
the skyline LU decomposition. The analysis does not renumber the equations, so that the LU without pivoting
is stable in the undrained steps. The VTK series add about 0.25 s per state (36 states in each of the two
runs, mostly the NeoPZ graph mesh of the integration point fields), about 18 s in all.

Files written in the working directory (the transposed-tangent run writes none):

| File | Contents |
|---|---|
| `embankment_history.csv` | t, λ, settlements at x = 0, 2, 4, 6 m, pp1, pp2 (Fig. 11); 11 undrained states (t = 0) and 25 consolidation steps |
| `embankment_convergence.csv` | normalized residual of every evaluation of every increment (stage 0 undrained, 1 consolidation; Sect. 6.6) |
| `embankment_nodal.scal_vec.{0,1,2}.vtk` | nodal displacement and pore pressure (native NeoPZ VTK) at the end of the loading, t = 10⁶ s and t = 10⁸ s (Fig. 12) |
| `embankment_gauss_{undrained,t1e6,t1e8}.vtk` | integration points: p', q, p_c, v0, type (0 elastic, 1 subcritical, 2 supercritical), σ' |
| `embankment_nodal_{undrained,t1e6,t1e8}.csv` | vertices: x, y, u_x, u_y, p, excess pore pressure p − γ_w(10 − y) (Fig. 12a, b, d) |
| `embankment_gauss_{undrained,t1e6,t1e8}.csv` | integration points: x, y, p', q, p_c, v0, type (Fig. 12c) |
| `embankment_elastic_history.csv`, `embankment_elastic_convergence.csv` | the same for the elastic variant |

The FLAC3D markers of Fig. 11 are the digitized histories of the Python package
(`dados/flac_historicos_digitalizados.json`); they are not written here.

## Viewing the solution in ParaView

The Cam-Clay and the elastic runs write all their 36 converged states with `mcc::TVTKSeries`, in
`vtk/embankment/` and `vtk/embankment_elastic/` (prefix `embankment` or `embankment_elastic`; 112 files and
about 18 MB in each directory):

| File | Contents |
|---|---|
| `embankment_nodal.vtk.series` → `embankment_nodal.scal_vec.<k>.vtk` | `Displacement` (vector) and `PorePressure` at the nodes, written by the NeoPZ graph mesh of the multiphysics mesh (`TPZAnalysis::DefineGraphMesh`/`PostProcess`) |
| `embankment_intpoints.vtk.series` → `embankment_intpoints.scal_vec.<k>.vtk` | the variables of the integration points of `TPZMatPoroElastoPlasticUP`, projected element by element on a discontinuous mesh by `TPZPostProcAnalysis` (as in the NeoPZ footing example): `MeanEffectiveStress` p', `DeviatoricStress` q, `PreconsolidationPressure` p_c, `PlasticType` (0 elastic, 1 subcritical, 2 supercritical), `VolumetricStrain` (tr ε, negative in compression), `SpecificVolume`, `EffectiveStressXX/YY/ZZ/XY`, `TotalStressXX/YY/ZZ/XY` (σ' − p_w I), `PrincipalEffectiveStress` (vector σ'1 ≥ σ'2 ≥ σ'3 at the points) and the tensors `EffectiveStress` and `TotalStress` |
| `embankment_gausspoints.vtk.series` → `embankment_gausspoints.<k>.vtk` | the 1800 integration points as a point cloud with their exact values: p', q, p_c, v0, type and σ' |
| `embankment_states.csv` | index k of each state with its time t (s) and load factor λ |

The time of the series is the index k of the state: 0 is the geostatic state, 1 to 10 are the undrained
increments (λ = 0.1 … 1, t = 0) and 11 to 35 the consolidation steps (t = 10^(2 + (k − 11)/4) s, up to
10⁸ s); `embankment_states.csv` gives t and λ of each k. Stresses are positive in tension (as in NeoPZ);
p' and q follow the soil mechanics convention (p' > 0 in compression). With 3 × 3 Gauss points the
projection is biquadratic: its values at the element vertices are the Lagrange extrapolation of the nine
Gauss values (the interpolation of `mcc::StressAtPoint`), and the jumps between neighbouring elements show
the discretization error. Being extrapolations, the vertex values can overshoot: the projected `PlasticType`
is not an integer (it ranges from about −1.9 to 4.4; the point cloud has the type of each point), q can be
slightly negative, and the projected principal stresses, sorted at the points, can be out of order where two
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
5. Excess pore pressure (Fig. 12a, b): select the reader (before the warp, so that the coordinates are the
   undeformed ones) and apply *Filters → Calculator* with *Result Array Name* `ExcessPorePressure` and the
   expression `PorePressure - 10*(10 - coordsY)`, i.e. p − γ_w (10 − y) kPa.
6. Integration point fields: open `vtk/embankment/embankment_intpoints.vtk.series` and color by
   `MeanEffectiveStress`, `DeviatoricStress`, `PreconsolidationPressure`, ... (the type of response of
   Fig. 12c is shown by the point cloud, item 7).
   To see them on the deformed mesh, apply *Filters → Resample With Dataset* (source: the nodal reader,
   destination: the integration point reader) to bring `Displacement` to its points, then *Warp By Vector*.
7. Point cloud: open `vtk/embankment/embankment_gausspoints.vtk.series`, set *Representation* to
   *Point Gaussian* (a *Gaussian Radius* of about 0.1 m) and color by `PlasticType`: at the last state the
   140 subcritical and 32 supercritical points of Fig. 12c appear.

The elastic run (`vtk/embankment_elastic/`) has the same files. The VTK files of the three states of Fig. 12
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
| `fig10_embankment_model` | Fig. 10: mesh, load, drained top, boundary conditions, settlement points and the zones pp1, pp2 | `embankment_nodal_undrained.csv` (vertices of the mesh) |
| `fig11_embankment_history` | Fig. 11: (a) settlements at x = 0, 2, 4, 6 m and (b) pore pressures pp1, pp2 from t = 0 to 10⁸ s, with the inset of pp2 against log t (Mandel–Cryer effect); markers: FLAC3D | `embankment_history.csv`, `reference/flac_historicos_digitalizados.json` |
| `fig12_embankment_fields` | Fig. 12: excess pore pressure at the end of the undrained loading (a) and at t = 10⁶ s (b), plastic integration points (c) and settlement (d) at t = 10⁸ s | `embankment_nodal_{undrained,t1e6,t1e8}.csv`, `embankment_gauss_t1e8.csv` |

The script draws the three states of the article; the fields of all 36 converged states are in the VTK series
(see *Viewing the solution in ParaView*).

## Results

Table 7 (settlements in m, positive downwards; pore pressures in kPa):

| | End of loading: this code | Python | article | FLAC3D | t = 10⁸ s: this code | Python | article | FLAC3D |
|---|---|---|---|---|---|---|---|---|
| Settlement, x = 0 | 0.152818 | 0.152818 | 0.153 | 0.140 | 0.275113 | 0.275113 | 0.275 | 0.193 |
| Settlement, x = 2 m | 0.152329 | 0.152329 | 0.152 | 0.135 | 0.269305 | 0.269305 | 0.269 | 0.186 |
| Settlement, x = 4 m | 0.067338 | 0.067338 | 0.067 | 0.055 | 0.164231 | 0.164231 | 0.164 | 0.104 |
| Settlement, x = 6 m | −0.038550 | −0.038550 | −0.039 | −0.042 | 0.022668 | 0.022668 | 0.023 | 0.004 |
| Pore pressure, pp1 | 33.0561 | 33.0561 | 33.1 | 18.1 | 5.0224 | 5.0224 | 5.0 | 5.1 |
| Pore pressure, pp2 | 55.7250 | 55.7250 | 55.7 | 62.4 | 25.1075 | 25.1075 | 25.1 | 25.1 |

Over the 36 monitored states the histories differ from the Python ones by at most 8·10⁻¹⁶ m (settlements)
and 8·10⁻¹⁴ kPa (pore pressures).

| Quantity | this code | Python / article |
|---|---|---|
| Nodes, pore pressure nodes | 661, 231 | article 661, 231 |
| Initial residual, largest component at the free equations | 8.8·10⁻¹³ kN/m | Python 8.4·10⁻¹³, article 8·10⁻¹³ |
| Initial vertical reaction of the base | 4600.000 kN/m | 23 × 20 × 10 = 4600 |
| p'0, q0 at the base | 84.0, 69.0 kPa | 84, 69 |
| Difference from FLAC3D at the end of the loading: settlement x = 0, heave x = 6 m | 9.2 %, 7.3 % (of the FLAC3D history value −0.0416 m) | article 9 %, 7 % |
| Reactions at t = 10⁸ s: x = 0, x = 20 m, base x, base y (kN/m) | 898.088, −795.210, −102.878, 4800.000 | Python 898.088, −795.210, −102.878, 4800.000; article 898.1, −795.2, −102.9, 4800.0 |
| Sum of the horizontal reactions | 9·10⁻⁷ kN/m | article < 10⁻⁵ |
| Integration points [elastic, sub, super], end of loading | [1800, 0, 0] | [1800, 0, 0] |
| Integration points [elastic, sub, super], t = 10⁸ s | [1628, 140, 32] | [1628, 140, 32] |
| Evaluations per increment, undrained | [5,5,5,4,4,4,4,4,4,4], mean 4.30 | identical (article 4.3) |
| Evaluations per increment, consolidation | [2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,4], mean 3.56 | identical (article 3.6) |
| Bisections | 0 | article: none |
| Residuals of the last undrained increment | 2.2e−2, 6.6e−4, 8.1e−7, 2.0e−12 | article: same |
| Residuals of the step ending at t = 10⁶ s | 7.7e−5, 5.6e−3, 1.7e−5, 2.0e−10 | article: same |
| Mandel–Cryer effect: largest pp2 | 58.840 kPa at t = 3.16·10⁵ s | Python 58.840 at 3.16·10⁵ s; article 55.7 → 58.8 kPa at 3.2·10⁵ s |
| pp2 at t = 10⁶ s | 55.815 kPa | Python 55.815; article: still at its undrained value |
| Settlement at x = 0, t = 10⁶ s, as a share of its consolidation part | 31.3 % | Python 31.3 %, article 31 % |
| Consolidation coefficient c_v = k(K + 4G/3) in zone pp2, end of loading and t = 10⁸ s | 1.24·10⁻⁶, 2.58·10⁻⁶ m²/s | Python data 1.24·10⁻⁶, 2.58·10⁻⁶; article 1–2.5·10⁻⁶ |
| c_v t/d² at t = 2.5·10⁵ s, d = 2.5 m | 0.050–0.103 | article 0.04–0.1 |
| Largest excess pore pressure, end of loading | 56.14 kPa at (1, 9) | Python 56.14 at (1, 9), article 56.1 |
| Largest excess pore pressure, t = 10⁶ s | 33.42 kPa at (0, 8) | Python 33.42 at (0, 8), article 33.4 |
| Heave of the top at t = 10⁸ s | 5.11 mm, from x = 9 m | Python 5.11 mm from x = 9 m; article up to 5 mm beyond x ≈ 9 m |
| Undrained loading with Dᵀ: evaluations per increment | [5,5,5,4,4,4,4,4,4,4]; residuals identical to those with D | Python `logUt`: identical |
| Elastic variant, final settlements x = 0, 2, 4, 6 m | 0.267524, 0.262540, 0.159354, 0.020157 | Python identical; article 0.268 at x = 0 |
| Plastic share of the final settlement at x = 0 | 0.0076 m | Python 0.0076, article 0.008 |

The c_v values are averages over the nine integration points of the pp2 element. The Python value is
computed in the same way from the states stored in `data_aterro.pkl`, with K = v0 p'/κ and
G = 3K(1 − 2ν)/(2(1 + ν)).

All the points stay elastic in the undrained loading, and the elastic tangent is symmetric because G is
evaluated at the start of the increment. The transposition therefore has no effect there, as in the
Python code.

The program prints these comparisons, with the reference values of the article and of the Python code
side by side.
