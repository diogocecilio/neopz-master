# SlopeRainfallA4

Caso A4 dos contraexemplos a Rahardjo et al. (2007) — "o lençol freático inicial é secundário" — com a **sucção
inicial sem limite** (o artigo limita a −75 kPa). Programa autocontido: todas as etapas estão explícitas em `main.cpp`
(malha geométrica do gmsh, malhas computacionais u P2 / p_w P1 / p_a P1 e multifísica, material
`TPZMatUnsaturatedUPwPa`, matriz estrutural e solver, análise `UnsatAnalysis`, contornos, estado inicial, chuva
unilateral, fator de segurança por redução de resistência e saída VTK).

Caso fixo (editar no topo de `main.cpp`): solo f50, α = 45°, H_s = 10 m, H_w = 10 m, chuva de 36 mm/h por 24 h,
48 h simuladas, hipóteses do artigo (percolação desacoplada, p_a = 0). Sucção inicial máxima 184 kPa; FS inicial
2,92 (2,40 com o limite de 75 kPa).

```
ninja SlopeRainfallA4
mkdir run_A4 && cd run_A4 && <build>/Projects/SlopeRainfallA4/SlopeRainfallA4
```
Saídas: `rain_A4.scal_vec.N.vtk` a cada hora (WaterPressure, AirPressure, Suction, Saturation, BishopPressure,
RelativePermeabilityWater, Displacement), `rain_A4.vtk.series` (abrir no ParaView: tempo em horas) e `fs_A4.csv`.
Tempo: ≈ 25–30 min (13 análises de FS).
