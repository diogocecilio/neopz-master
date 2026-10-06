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
./AbaqusTriaxialConsolidation            # all parts, about 1.5 min
./AbaqusTriaxialConsolidation mp axi     # material point and axisymmetric models, about 5 s
```

Files: `abaqus_<run>.csv` (δ/H, p'_A, q_A, σ_a, max p_w, ε_v), `abaqus_<run>.vtk` and `abaqus_<run>_gauss.vtk`
(final nodal fields and integration points), `abaqus_material_point_*.csv`, `abaqus_closed_form_*.csv` and
`abaqus_table8.csv`.

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
| Largest excess pore pressure (kPa) | 3.8·10⁻⁴ | 4·10⁻⁴ |
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
