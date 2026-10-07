# Examples of the Modified Cam-Clay u–p article

This folder contains one project per numerical example of the article

> D. Lira Cecílio, *Return mapping for Modified Cam-Clay plasticity in rotated Haigh–Westergaard space
> with consistent tangent operator and coupled u–p consolidation*.

The examples are the C++/NeoPZ counterparts of the Python transcription of the Wolfram Language packages
(`camclay_hw.py` = `camclay-perf.m`, `fe_user.py` = `poro-camclay-fem.m`, drivers `gen_data.py`,
`gen_data3d.py`, `aterro_elastic.py` and `fig_surface.py`) and of the rival integration schemes that existed only
in Python (`rivais.py`: implicit backward-Euler variants and adaptive explicit Runge-Kutta schemes, ported in
IntegrationSchemes). They reproduce the numbers of the Python code to round-off and the tables of the article.

| Project | Article | Python driver |
|---|---|---|
| [YieldSurfaceProjection](YieldSurfaceProjection) | Figs. 1 and 2 | `fig_surface.py` |
| [TaylorTest](TaylorTest) | Sect. 4.5, Fig. 3; Table 10 (Taylor slopes of the tangent operators) | `gen_data.py taylor`, `tangentes` (Taylor slopes) |
| [RS2Triaxial](RS2Triaxial) | Sect. 6.1, Figs. 4 and 5, Table 2 | `gen_data.py rs2` |
| [FrozenBulkModulus](FrozenBulkModulus) | Sect. 6.1, Table 3 | `gen_data.py frozen` |
| [IntegrationSchemes](IntegrationSchemes) | Sect. 6.2, Fig. 6, Table 4 | `gen_data.py rivais`, `rivais.py` |
| [FLAC3DTriaxial](FLAC3DTriaxial) | Sects. 6.3 and 6.7, Figs. 4 and 7, Tables 5 and 10 | `gen_data.py itasca` |
| [TerzaghiConsolidation](TerzaghiConsolidation) | Sect. 6.4, Fig. 8, Table 6 | `gen_data.py terzaghi`, `gen_data3d.py` |
| [AbaqusTriaxialConsolidation](AbaqusTriaxialConsolidation) | Sects. 6.5 and 6.7, Figs. 9–11, Tables 7, 9 and 10 | `gen_data.py abaqus abaqus_mp abaqus_states`, `gen_data3d.py` |
| [EmbankmentConsolidation](EmbankmentConsolidation) | Sects. 6.6 and 6.7, Figs. 12–14, Tables 8 and 10 | `gen_data.py aterro`, `aterro_elastic.py` |

Section, table and figure numbers of the article v0.7 (the figure of the single-element model, Fig. 4, was added,
so the later figures moved by one with respect to v0.6). The material-point examples
(YieldSurfaceProjection, TaylorTest, FrozenBulkModulus, IntegrationSchemes and the material-point tests of
RS2Triaxial) have no finite element mesh.

## Library classes (Material/Plasticity)

| Class | Role | Article / WL routines |
|---|---|---|
| `TPZYCModifiedCamClayRHW` | yield function, hardening law, local residuals, Jacobian, Newton projection and Jacobian of the projection | (11)–(13), (17)–(22), (A.1); `HardeningCC`, `PhiCC`, `ResCC`, `JacCC`, `dResdTrialCC`, `ProjectHWCC`, `GradCC` |
| `TPZPlasticStepModifiedCamClay` | elastic predictor with the porous law, spectral decomposition, stress update and consistent tangent (also a linear elastic option); tangent returned to the global iterations selected by `SetTangentMode`: consistent `D` (default), `D^T`, `(D+D^T)/2`, continuum operator, central differences (Sect. 6.7, Table 10) | Algorithm 1, (8)–(10), (14)–(16); `TrialStressCC`, `ProjectStressCC`, `ComputedDep`; `TANGENT['mode']` of `gen_data.py`, `continuum_tangent` |
| `TPZMatPoroElastoPlasticUP` | multiphysics u–p material with memory, plane strain / axisymmetry / 3D with the six-row operator; post-processing variables of the integration points (p', q, p_c, type of response, stresses, ...) read from the memory | (24)–(28); `ComputeBN`, `ContributePorous`, `ContributePlasticity` |
| `TPZPoroElastoPlasticUPAnalysis` | incremental Newton driver with time step, load factor, controlled displacement, bisection and reactions; convergence records of the converged increments (`StepLog`) and work counters with the failed attempts (`NGlobalIterations`, `NBisections`) | Sect. 5.4; `SolveStepUP`, `AdvanceUP`, `IterativeProcessUP`, `ReactionByMarker`; `itcount`, `ncut` of `fe_user.py` |

`Common/MCCPaperTools.h` gathers the utilities shared by the projects: structured Q8–Q4 and Hex20–Hex8 meshes
(box, block, unit cube, slab and the quarter cylinder with quadratic geometry), CSV files of the geometric mesh
for plotting it in Python (`mcc::WriteMeshCSV`: nodes, elements, boundary faces and edges with the mid-edge nodes
of the curved geometry), the displacement (serendipity), pore pressure and
multiphysics meshes, output at the integration points, interpolation of the stress at a point, the VTK file
series of every converged state for ParaView (`mcc::TVTKSeries`), material point drivers and the closed-form
solutions of Appendix B.

Each project is written as a class in the style of the NeoPZ examples (the `Footing` class): geometric mesh,
computational meshes, analysis with structural matrix (`TPZSkylineNSymStructMatrix`) and solver
(`TPZStepSolver`, LU), incremental solution and post-processing (CSV files with the data of the figures and VTK
files of the nodal fields and of the integration points).

## Building and running

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DCMAKE_BUILD_TYPE=Release <path to neopz>
make -j 4
./Projects/AbaqusTriaxialConsolidation/AbaqusTriaxialConsolidation
```

The executables write their files in the current directory and print the comparison of their results with
the values of the article and of the Python code. Run times (Release build, one core): the material-point examples
take less than half a second (YieldSurfaceProjection 0.4 s, TaylorTest 0.07 s, FrozenBulkModulus 0.4 s,
IntegrationSchemes 0.3 s); the finite element examples, all three-dimensional with Hex20–Hex8 elements, take from
a few seconds (TerzaghiConsolidation, FLAC3DTriaxial, RS2Triaxial) to a few minutes (AbaqusTriaxialConsolidation,
EmbankmentConsolidation, whose comparison of the tangent operators of Table 10 repeats the analysis with each
operator). The README of each project gives its run times and the arguments that select its parts (for example
`novtk`, which skips the VTK series; see below). The documentation of the classes is generated with
`-DBUILD_DOCS=ON` (Doxygen group *Examples of the Modified Cam-Clay u-p article*).

## Figures

Each project has a script `plot_figures.py` (Python 3 with numpy and matplotlib) that draws the figures of the
article from the CSV files of its executable, with the style of the figures of the article
(`Common/mcc_figstyle.py`) and the digitized reference curves of the Python code (`<project>/reference/`):

```
cd <run directory>          # where the executable wrote its CSV files
./FLAC3DTriaxial
python3 <neopz>/Projects/FLAC3DTriaxial/plot_figures.py      # writes figures/fig07_flac3d_triaxial.pdf/.png
```

| Project | Files (PDF and PNG in `<run directory>/figures`) |
|---|---|
| YieldSurfaceProjection | `fig01_mcc_surface`, `fig02_meridian_projection` |
| TaylorTest | `fig03_taylor_test` (and the same figure for the `std::mt19937_64` sample); `supplementary_taylor_operators` (Taylor slopes of Table 10; not a figure of the article) |
| RS2Triaxial | `fig05_rs2_triaxial`; `supplementary_rs2_element_check` (not a figure of the article) |
| FrozenBulkModulus | `supplementary_table03_frozen_bulk_modulus` (Table 3; not a figure of the article) |
| IntegrationSchemes | `fig06_integration_schemes`; `supplementary_secant_first_increment` (single secant step in the first increment of the drained test with OCR = 10; not a figure of the article) |
| FLAC3DTriaxial | `fig04_single_element_model` (the Hex20–Hex8 element of the RS2 and FLAC3D tests with its boundary conditions), `fig07_flac3d_triaxial` |
| TerzaghiConsolidation | `fig08_terzaghi_consolidation` |
| AbaqusTriaxialConsolidation | `fig09_abaqus_model`, `fig10_abaqus_states`, `fig11_abaqus_results`; `supplementary_abaqus_softening`, `supplementary_abaqus_tangents` (not figures of the article) |
| EmbankmentConsolidation | `fig12_embankment_model`, `fig13_embankment_history`, `fig14_embankment_fields` |

The numbers in the file names are the figure numbers of the article v0.7.

The scripts are ports of `fig_surface.py` and `figs.py` of the Python code, reading the CSV files instead of the
pickled results; the figures are the same as those of the article. The font of the article (TeX Gyre Heros) is
used when it is installed; otherwise a font with the metrics of Helvetica (Nimbus Sans, Liberation Sans or
FreeSans) or DejaVu Sans.

## Viewing the solution in ParaView

EmbankmentConsolidation and AbaqusTriaxialConsolidation write the solution of every converged state (every
increment or time step) as VTK file series, one directory per run under `vtk/` (`mcc::TVTKSeries`):

* `<prefix>_nodal.vtk.series`: displacement and pore pressure at the nodes (NeoPZ graph mesh of the
  multiphysics mesh);
* `<prefix>_intpoints.vtk.series`: the variables of the integration points (p', q, p_c, type of response,
  volumetric strain, effective and total stresses, principal stresses), projected element by element on a
  discontinuous mesh by `TPZPostProcAnalysis`, as in the NeoPZ footing example;
* `<prefix>_gausspoints.vtk.series`: the integration points as a point cloud;
* `<prefix>_states.csv`: the time, load factor or displacement of each state.

Open the `.vtk.series` files in ParaView (5.5 or later): the time controls then show the time of each state
(the state index for the embankment, δ/H for the triaxial test). The READMEs of the two projects describe the
files and the steps in ParaView (field, time, Warp By Vector, excess pore pressure with the Calculator,
point cloud). The argument `novtk` skips these files; the READMEs of the two projects give the run times and
the size of the files with them.
