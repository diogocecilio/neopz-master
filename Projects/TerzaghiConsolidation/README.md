# TerzaghiConsolidation — Sect. 6.3, Fig. 6 and Table 5

Consolidation of a linear elastic column (H = 10 m, E = 10⁴ kPa, ν = 0.25, k = 10⁻⁶ m²/(kPa·s),
c_v = 0.012 m²/s) loaded with q = 10 kPa on the drained top; base fixed and impermeable, lateral faces
impermeable and restrained horizontally; incompressible constituents (α_B = 1, 1/M_B = 0).

* 2D: 1 × 10 Q8–Q4 elements in plane strain, 3 × 3 Gauss points.
* 3D: 1 × 1 × 10 Hex20–Hex8 elements, 3 × 3 × 3 Gauss points (the article reports identical results).
* Load applied in an undrained step (Δt = 0) followed by 101 time steps (20 per decade) from
  T = c_v t/H² = 10⁻⁵ to 1, plus T = 0.001, 0.01, 0.1 and 0.5. The elastic response is given by
  `TPZPlasticStepModifiedCamClay::SetLinearElastic`.

## Running

```
./TerzaghiConsolidation
```

Run time: about 3 s. Files: `terzaghi_q8q4_history.csv` and `terzaghi_hex20hex8_history.csv` (time,
settlement of the top and pore pressures at the vertices of x = 0), `terzaghi_fig6a.csv` (p_w/q along the
column for the four times, numerical and exact), `terzaghi_fig6b.csv` (degree of consolidation),
`terzaghi_table5.csv`, and VTK files of the nodal fields at the four times.

## Results (Table 5)

| T | 0.001 | 0.01 | 0.1 | 0.5 |
|---|---|---|---|---|
| max \|p_w/q − exact\| (this code / article) | 0.075 / 0.075 | 0.011 / 0.011 | 0.007 / 0.007 | 0.015 / 0.015 |
| settlement, this code (mm) | 0.366 | 0.954 | 2.959 | 6.289 |
| settlement, article (mm) | 0.366 | 0.954 | 2.959 | 6.289 |
| settlement, exact (mm) | 0.297 | 0.940 | 2.974 | 6.366 |

Undrained step: p_w/q = 1.27 at y = 9 m and 0.93 at y = 8 m (article: 1.27 and 0.93). The histories differ
from those of the Python code by less than 6·10⁻¹⁵ m in the settlement and 5·10⁻¹¹ kPa in the pore
pressures, and the Hex20–Hex8 model differs from the Q8–Q4 model by 2.5·10⁻¹⁴ (article: 5·10⁻¹¹).
