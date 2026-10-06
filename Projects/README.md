# Examples of the Modified Cam-Clay u–p article

This folder contains one project per numerical example of the article

> D. Lira Cecílio, *Return mapping for Modified Cam-Clay plasticity in rotated Haigh–Westergaard space
> with consistent tangent operator and coupled u–p consolidation*.

The examples are the C++/NeoPZ counterparts of the Python transcription of the Wolfram Language packages
(`camclay_hw.py` = `camclay-perf.m`, `fe_user.py` = `poro-camclay-fem.m`, drivers `gen_data.py`,
`gen_data3d.py`, `aterro_elastic.py` and `fig_surface.py`). They reproduce the numbers of the Python code to
round-off and the tables of the article.

| Project | Article | Python driver |
|---|---|---|
| [YieldSurfaceProjection](YieldSurfaceProjection) | Figs. 1 and 2 | `fig_surface.py` |
| [TaylorTest](TaylorTest) | Sect. 4.5, Fig. 3 | `gen_data.py taylor` |
| [RS2Triaxial](RS2Triaxial) | Sect. 6.1, Fig. 4, Table 2 | `gen_data.py rs2` |
| [FrozenBulkModulus](FrozenBulkModulus) | Sect. 6.1, Table 3 | `gen_data.py frozen` |
| [FLAC3DTriaxial](FLAC3DTriaxial) | Sect. 6.2, Fig. 5, Table 4 | `gen_data.py itasca` |
| [TerzaghiConsolidation](TerzaghiConsolidation) | Sect. 6.3, Fig. 6, Table 5 | `gen_data.py terzaghi`, `gen_data3d.py` |
| [AbaqusTriaxialConsolidation](AbaqusTriaxialConsolidation) | Sects. 6.4 and 6.6, Figs. 7–9, Tables 6 and 8 | `gen_data.py abaqus abaqus_mp abaqus_states`, `gen_data3d.py` |
| [EmbankmentConsolidation](EmbankmentConsolidation) | Sect. 6.5, Figs. 10–12, Table 7 | `gen_data.py aterro`, `aterro_elastic.py` |

## Library classes (Material/Plasticity)

| Class | Role | Article / WL routines |
|---|---|---|
| `TPZYCModifiedCamClayRHW` | yield function, hardening law, local residuals, Jacobian, Newton projection and Jacobian of the projection | (11)–(13), (17)–(22), (A.1); `HardeningCC`, `PhiCC`, `ResCC`, `JacCC`, `dResdTrialCC`, `ProjectHWCC`, `GradCC` |
| `TPZPlasticStepModifiedCamClay` | elastic predictor with the porous law, spectral decomposition, stress update and consistent tangent (also a linear elastic option) | Algorithm 1, (8)–(10), (14)–(16); `TrialStressCC`, `ProjectStressCC`, `ComputedDep` |
| `TPZMatPoroElastoPlasticUP` | multiphysics u–p material with memory, plane strain / axisymmetry / 3D with the six-row operator | (24)–(28); `ComputeBN`, `ContributePorous`, `ContributePlasticity` |
| `TPZPoroElastoPlasticUPAnalysis` | incremental Newton driver with time step, load factor, controlled displacement, bisection and reactions | Sect. 5.4; `SolveStepUP`, `AdvanceUP`, `IterativeProcessUP`, `ReactionByMarker` |

`Common/MCCPaperTools.h` gathers the utilities shared by the projects: structured Q8–Q4 and Hex20–Hex8 meshes
(including the quarter cylinder with quadratic geometry), the displacement (serendipity), pore pressure and
multiphysics meshes, output at the integration points, interpolation of the stress at a point, material point
drivers and the closed-form solutions of Appendix B.

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
the values of the article and of the Python code. Run times (Release build, one core): YieldSurfaceProjection
0.3 s, TaylorTest 0.05 s, RS2Triaxial 1.9 s, FrozenBulkModulus 0.4 s, FLAC3DTriaxial 1.6 s,
TerzaghiConsolidation 2.9 s, AbaqusTriaxialConsolidation 68 s and EmbankmentConsolidation 35 s. The documentation of the classes is generated with `-DBUILD_DOCS=ON` (Doxygen
group *Examples of the Modified Cam-Clay u-p article*).
