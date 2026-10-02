# Adensamento u-p 3D com Cam-Clay modificado no NeoPZ

Porte para o NeoPZ do FE u-p (Biot) 3D do Python (`python/fe3d_up.py`) com o Cam-Clay modificado
`TPZModifiedCamClay` (`Material/Plasticity`), e dois exemplos de verificação:

| executável | problema | referência Python |
|---|---|---|
| `AterroCamClay` | aterro sobre fundação Cam-Clay (exemplo do FLAC3D, Itasca) | `aterro_itasca.py`, `aterro_elementos.py` |
| `TriaxialAbaqusCamClay` | Abaqus Benchmarks 1.15.2, adensamento de um corpo de prova triaxial, hexaedros de 20 nós | `triaxial_abaqus.py` |

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON <neopz>
make AterroCamClay TriaxialAbaqusCamClay
cd Projects2/PoroCamClay3D
./AterroCamClay                  # hex20 e hex8 (ou: ./AterroCamClay hex8)
./TriaxialAbaqusCamClay          # placas lisa e rugosa, hex20, malha 2x2x4, 150 passos
./TriaxialAbaqusCamClay rugosa hex20 3 3 8 150       # placa, elemento, nc nr nz, passos
python3 <neopz>/Projects2/PoroCamClay3D/compara_python.py   # opcional: compara com o Python
```

## Estrutura (nativa do NeoPZ)

Os dois exemplos usam as estruturas do NeoPZ, como `Projects2/PoroElastic` e `Projects2/GeoMecDeterm`:

* **Material** `TPZMatPoroElastoPlastic3DMem<T, TMEM>` (`Material/Plasticity`, instanciado para
  `TPZModifiedCamClay`): material multifísico `TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
  TPZMatWithMem<TPZElastoPlasticMem>>`, versão 3D com Newton do `TPZMatPoroElastoPlastic2DMem`.
  `datavec[0]` = u (H1 vetorial), `datavec[1]` = p (H1 escalar). Formulação total, Euler implícito:

  ```
  R_u = ∫ Bᵀσ'(ε) - α ∫ Bᵀm N_p p - ∫ N_uᵀ b - ∫ N_uᵀ t
  R_p = ∫ N_pᵀ [α (ε_v - ε_v,n) + S (p - p_n)] + Δt ∫ ∇N_pᵀ (k/μ)(∇p - ρ_f g)
  ek  = [[∫ Bᵀ D_ep B, -Q], [Qᵀ, S M_p + Δt H]],   ef = -R
  ```

  σ' vem de `TPZModifiedCamClay::ApplyStrainComputeSigma` com o estado da memória (ε_n, ε^p_n, α_n e σ'_n,
  necessária nas leis hipoelásticas). Condições de contorno (Val2 = {v_x, v_y, v_z, p}): 0 u = v,
  1 tração, 2 p, 3 u_i = 0 nas direções marcadas, 5 pressão normal, 6 u_i prescrito nas direções de Val1,
  12 tração + p, 16 tipo 6 + p (Dirichlet por penalidade). Variáveis: `Displacement`, `Pressure`,
  `ExcessPressure`, `Flux`, `PressureGradient`. Modos de montagem: completo, só forças externas (norma de
  referência do critério de convergência) e sem penalidade (reações, p.ex. a força na placa do triaxial).
* **Malhas** (`PoroCamClayNativo.cpp`):
  * `cmeshU`: H1 vetorial com memória (`SetAllCreateFunctionsContinuousWithMem`), material
    `TPZMatCamClayPostProc` (um `TPZMatElastoPlastic<TPZModifiedCamClay>`);
  * `cmeshP`: H1 escalar de ordem 1 (`TPZNullMaterial`);
  * `TPZMultiphysicsCompMesh` com `BuildMultiphysicsSpace({1, 1}, {cmeshU, cmeshP})`.

  O material u-p usa os índices de memória do elemento atômico de u (`datavec[0].intGlobPtIndex`), e as duas
  malhas compartilham o mesmo vetor de memória (`fMat->GetMemory() = fMatU->GetMemory()`): o que o material u-p
  grava nos pontos de integração (σ', ε^p, α, p) é lido pelo material de u no pós-processamento.
* **Análise**: `TPZLinearAnalysis(mphys, true)` (renumeração de banda) com
  `TPZSkylineNSymStructMatrix` + `ELU` (o sistema u-p é não simétrico). Cada passo é um laço de Newton com
  `Assemble()` / `Solve()` / `LoadSolution()`; ao convergir, `SetUpdateMem(true)` + `AssembleResidual()`
  grava o estado na memória. Em caso de falha (Newton ou return mapping) o passo é subdividido.
* **Saída VTK**:

  ```cpp
  an.DefineGraphMesh(3, {"Pressure", "ExcessPressure"}, {"Displacement", "Flux"}, base + "_up.vtk");
  an.SetStep(n);  an.PostProcess(res);                       // <base>_up.scal_vec.n.vtk

  TPZPostProcAnalysis post;  post.SetCompMesh(cmeshU);       // variáveis dos pontos de integração
  post.SetPostProcessVariables(matids, vars);
  post.DefineGraphMesh(3, escalares, vetores, tensores, base + "_tensoes.vtk");
  post.TransferSolution();  post.SetStep(n);  post.PostProcess(res);   // <base>_tensoes.scal_vec.n.vtk

  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, arquivo, true);       // <base>_malha.vtk
  ```

  `_up`: `Pressure`, `ExcessPressure`, `Displacement`, `Flux`. `_tensoes`: `PorePressure`,
  `ExcessPorePressure`, `MeanEffectiveStress` (p'), `DeviatoricStress` (q), `VolHardening` (α),
  `PreconsolidationPressure` (p_c), `PlasticPoint`, `Displacement`, `EffectiveStress` e `TotalStress`
  (tensores). Os campos de `_tensoes` são a projeção L2 dos valores nos pontos de integração feita pelo
  `TPZPostProcAnalysis`, então podem ultrapassar um pouco os valores extremos nos cantos (p.ex. q < 0 junto
  à quina da placa rugosa). Os arquivos `*.scal_vec.N.vtk` abrem no ParaView como série temporal.

### Elementos

* `hex8`: `TPZGeoCube`, u e p lineares (Q1-Q1, 2x2x2) — o mesmo elemento do Python.
* `hex20`: `TPZQuadraticCube` (20 nós, nós de meio de aresta sobre o arco no triaxial), u com o H1 de ordem 2
  do NeoPZ e p linear, integração 3x3x3. No hexaedro o H1 de ordem 2 do NeoPZ é o Q2 hierárquico
  (27 funções: vértices, arestas, faces e bolha), não o serendipity de 20 funções do Python/Abaqus; a
  diferença de resultado é pequena (ver abaixo). Com Q2 a integração reduzida 2x2x2 (C3D20RP/CAX8RP do
  Abaqus) deixa modos espúrios de energia nula e a matriz fica singular, por isso ela não é oferecida aqui.

### Versão de referência com o elemento serendipity (`serendipity/`)

A primeira versão do porte montava o sistema fora do `TPZCompMesh` (`serendipity/TPZPoroCamClayUP`), com as
funções do próprio mapeamento `TPZQuadraticCube` (serendipity de 20 nós), as regras `TPZIntCube3D` e um
escritor VTK próprio, para reproduzir o Python e o Abaqus elemento por elemento (inclusive o `hex20r`). Ela
continua disponível como referência: `AterroCamClaySerendipity` e `TriaxialAbaqusCamClaySerendipity`
(saídas com o prefixo `serendipity_`).

## Resultados

**Aterro (FLAC3D)**, malha 20x1x10, 10 incrementos de carga não drenados + 25 passos de adensamento até
10⁸ s; saída `aterro_<elem>_up/_tensoes.scal_vec.0..35.vtk` e `aterro_<elem>_historico.csv`:

* hex8: NeoPZ e Python coincidem em todos os instantes (diferenças < 5e-11 m nos recalques e < 5e-9 kPa nas
  poropressões, no limite da precisão dos CSV), com o equilíbrio inicial exato (|R_u| livre = 1.8e-13 kN);
* hex20 (u Q2) x hex20 serendipity do Python: recalque em x = 0 igual em 4 dígitos (0.1528 m no fim da etapa
  não drenada e 0.2751 m em 10⁸ s), poropressões pp1/pp2 iguais a 0.01 kPa; a maior diferença é no recalque
  em x = 4 m, na borda da carga (0.1674 x 0.1642 m em 10⁸ s);
* FLAC3D: ≈ 0.14 m no fim da etapa não drenada e ≈ 0.19 m em 10⁸ s. A diferença do recalque de adensamento
  já está no modelo Python (é dominada pela lei elástica, ver o estudo de sensibilidade em
  `aterro_itasca.py`), não vem do porte.

**Abaqus 1.15.2**, quarto de cilindro 2x2x4 hex20 (321 nós, 1742 equações), 150 passos até δ/H = 0.6 em
400 dias; saída `triaxial_<placa>_hex20_224_up/_tensoes.scal_vec.0..15.vtk` e `.csv` (δ/H, p' e q no ponto A,
σ_a sob a placa, |p| máx):

| q no ponto A (kPa), δ/H = 0.51 / 0.60 | placa lisa | placa rugosa |
|---|---|---|
| NeoPZ hex20 (u Q2) | 148.8 / 149.5 | 145.8 / 142.7 |
| Python hex20 serendipity | 148.8 / 149.5 | 145.2 / 141.8 |
| Python hex20r serendipity (integração reduzida) | 148.8 / 149.5 | 154.5 / 155.9 |
| Abaqus CAX8RP (digitalizado, δ/H = 0.51) | ≈ 148 | ≈ 153 |

* placa lisa: estado homogêneo; NeoPZ (Q2) e Python (serendipity) diferem < 2e-4 kPa em todos os passos e
  ficam a 0.1 kPa da solução homogênea (600 incrementos);
* placa rugosa: o Q2 com integração completa fica com o serendipity completo (diferença ≤ 1 kPa); o
  amolecimento depois do pico vem da integração completa junto à quina da placa — o Abaqus usa integração
  reduzida, e o hex20r serendipity (`TriaxialAbaqusCamClaySerendipity`) reproduz a curva do Abaqus.

`compara_python.py` refaz os casos com o código Python de `python/`, compara todos os passos das versões
nativa e serendipity e gera `comparacao_aterro.png` e `comparacao_triaxial.png`.
