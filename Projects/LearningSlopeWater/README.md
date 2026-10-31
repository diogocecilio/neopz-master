# LearningSlopeWater: talude saturado, transiente e FS ao longo do tempo

Etapa "com água" do estudo do Mohr–Coulomb PV2 com o material não saturado u–p_w–p_a (`TPZMatUnsaturatedUPwPa`),
ar passivo (p_a = 0) e limite saturado do van Genuchten (a_VG = 10⁻⁸: S_w = 1, k_rw = 1).

* Talude de de Souza Neto et al. (Fig. 8.30): h = 10 m, 45°, E = 20 MPa, ν = 0,49, c = 50 kPa, φ = ψ = 20°,
  γ = 20 kN/m³; k_s = 10⁻⁶ m/s. Malha do Gmsh de `LearningMohrCoulomb/meshes/talude.msh` (P2/P1/P1).
* Estado inicial: hidrostático com o nível d'água no pé (y = 30 m) e estado geostático com essas pressões.
* Transiente: em t = 0 o nível no contorno esquerdo sobe para 38 m; no direito fica em 30 m; o terreno a jusante
  do pé fica drenado (p_w = 0); base, crista e face impermeáveis. 50 passos em progressão geométrica até 2·10⁶ s.
* FS por redução de resistência (SRM) com as pressões do instante congeladas, a cada 10 passos; a memória e a
  solução do transiente são restauradas depois de cada SRM.
* Dois casos: `acoplado` (deformação no balanço de massa da água) e `desacoplado`
  (`SetHydroMechanicalCoupling(false)`: só a compressibilidade da água armazena).

Saídas: `fs_tempo.csv` (caso, passo, t, FS) e `pressao_tempo.csv` (caso, passo, t, ponto, x, y, p_w).
Gráficos: `python3 plot_fs_pressao.py fs_tempo.csv pressao_tempo.csv fs_pressao.png`.

## Resultados (31/10/2026)

| t (h) | 0 | 0,41 | 2,8 | 16,6 | 96 | 556 |
|---|---|---|---|---|---|---|
| FS acoplado | 2,481 | 2,481 | 2,359 | 2,209 | 2,200 | 2,200 |
| FS desacoplado | 2,481 | 2,266 | 2,200 | 2,200 | 2,200 | 2,200 |

* O regime permanente é o mesmo nos dois casos (p_w idêntico em todos os pontos de controle, FS = 2,200).
* Sem acoplamento o armazenamento é só n/K_w e o permanente chega em ~3 h; com acoplamento entra a
  compressibilidade do esqueleto e o processo segue a escala do adensamento (~20 h).
* FS inicial 2,48 acima do talude seco (2,275): no limite saturado a pressão negativa acima do nível d'água
  entra com χ = 1 e aumenta a tensão efetiva; a subida do nível reduz o FS a 2,20.
