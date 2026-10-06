# FLAC3DTriaxial — Sect. 6.2, Fig. 5 and Table 4

Drained and undrained triaxial tests of the FLAC3D verification problem on a single axisymmetric
Q8–Q4 u–p element (1 m × 1 m, 2 × 2 Gauss points), from an isotropic effective stress p'0 = 5 kPa with
overconsolidation ratios R = p'c0/p'0 = 1.6 (subcritical) and R = 8 (supercritical).

* Material (Table 1): M = 1.02, λ = 0.2, κ = 0.05, v_λ = 3.32 (v0 on the NCL), G = 250 kPa, porous law.
* Boundary conditions: u_r = 0 on the axis, u_z = 0 at the base, total cell pressure on the lateral face,
  vertical displacement of the top controlled.
* Drained tests: p_w = 0 at the four vertices, ε_a up to 50 % in 500 increments.
* Undrained tests: no flow (k = 0, Δt = 0), pore fluid with K_w = 2·10⁴ kPa (M_B = K_w/n), ε_a up to 10 %
  in 400 increments.
* Sect. 6.6: drained test with R = 1.6 in 50 increments with the consistent operator D and with its
  transpose Dᵀ (the transposition has no effect in this homogeneous test).

The class `FLAC3DTriaxial` follows the structure of the NeoPZ examples: `CreateGeoMesh`,
`CreateCompMesh` (displacement, pore pressure and multiphysics meshes with the material
`TPZMatPoroElastoPlasticUP` and its boundary conditions), `Run` (analysis `TPZPoroElastoPlasticUPAnalysis`
with `TPZSkylineNSymStructMatrix` and `TPZStepSolver` LU, increments and post-processing).

## Running

```
./FLAC3DTriaxial
```

Run time: about 2 s. Files: `flac3d_<test>.csv` (ε_a, p', q, v, u for Fig. 5), `flac3d_<test>.scal_vec.0.vtk`
(nodal fields) and `flac3d_<test>_gauss.vtk` (state of the integration points).

## Figures

```
python3 <neopz>/Projects/FLAC3DTriaxial/plot_figures.py [run directory] [-o output directory]
```

Run it after the executable, with the directory of its CSV files (default: the current directory). The figures
are written as PDF and PNG to `<run directory>/figures` (or to the output directory); Python 3 with numpy and
matplotlib is needed.

| File | Article | Data |
|---|---|---|
| `fig05_flac3d_triaxial` | Fig. 5: drained (a, b) and undrained (c, d) tests with R = 1.6 and R = 8, q–ε_a and stress paths p'–q with the CSL and the initial yield surface; dashed lines: closed-form solutions; squares: FLAC3D final states (Table 4) | `flac3d_<drained\|undrained>_R<1.6\|8>.csv` |

The script also prints the final states, the peaks and the closed-form values at the same ε_a.

## Results (final states, p', q and u in kPa)

| Test | | this code | Python | article (this work) | FLAC3D |
|---|---|---|---|---|---|
| Drained, R = 1.6 | p' / q / v | 7.57329 / 7.71988 / 2.81120 | 7.57329 / 7.71988 / 2.81120 | 7.573 / 7.720 / 2.811 | 7.573 / 7.718 / 2.811 |
| Drained, R = 8 | p' / q / v | 7.58388 / 7.75165 / 2.81051 | 7.58388 / 7.75165 / 2.81051 | 7.584 / 7.752 / 2.811 | 7.583 / 7.747 / 2.811 |
| Undrained, R = 1.6 | p' / q / u | 4.23438 / 4.31853 / 2.20513 | 4.23438 / 4.31853 / 2.20513 | 4.234 / 4.319 / 2.205 | 4.234 / 4.312 / 2.203 |
| Undrained, R = 8 | p' / q / u | 14.0479 / 14.4224 / −4.24044 | 14.0479 / 14.4224 / −4.24044 | 14.048 / 14.422 / −4.240 | 14.05 / 14.42 / −4.241 |

Peaks with R = 8: q = 18.108 kPa at ε_a = 3.0 % (drained) and 15.048 kPa at ε_a = 4.0 % (undrained), as in the
article. Evaluations of the residual per increment: 2.228 and 2.34 (drained), 2.0025 (undrained), identical
to the Python code; with D and Dᵀ (50 increments): 3.02 and 3.02.
