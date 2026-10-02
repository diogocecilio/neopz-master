# Adensamento u-p 3D com Cam-Clay modificado no NeoPZ

Porte para C++/NeoPZ do FE u-p (Biot) 3D do Python (`python/fe3d_up.py`) com o Cam-Clay modificado
`TPZModifiedCamClay` (`Material/Plasticity`), e dois exemplos de verificação:

| executável | problema | referência Python |
|---|---|---|
| `AterroCamClay` | aterro sobre fundação Cam-Clay (exemplo do FLAC3D, Itasca) | `aterro_itasca.py`, `aterro_elementos.py` |
| `TriaxialAbaqusCamClay` | Abaqus Benchmarks 1.15.2, adensamento de um corpo de prova triaxial, hexaedros de 20 nós | `triaxial_abaqus.py` |

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON <neopz>
make AterroCamClay TriaxialAbaqusCamClay
cd Projects2/PoroCamClay3D
./AterroCamClay                 # hex20 e hex8 (ou: ./AterroCamClay hex20)
./TriaxialAbaqusCamClay         # placas lisa e rugosa, hex20 e hex20r, malha 2x2x4, 150 passos
./TriaxialAbaqusCamClay rugosa hex20 3 3 8 150      # placa, elemento, nc nr nz, passos
python3 <neopz>/Projects2/PoroCamClay3D/compara_python.py   # opcional: compara com o Python
```

## Solver (`TPZPoroCamClayUP`)

* Malha: `TPZGeoMesh` com `TPZGeoCube` (8 nós) ou `TPZQuadraticCube` (20 nós) e faces de contorno
  `TPZGeoQuad` / `TPZQuadraticQuad` cujo material id é o marcador da face (`malhas.cpp`: caixa estruturada e
  quarto de cilindro "O-grid" com os nós de meio de aresta no arco).
* Elementos: `hex8` (Q1-Q1, 2x2x2), `hex20` (u serendipity quadrático, p trilinear nos vértices: Q2-Q1,
  3x3x3) e `hex20r` (hex20 com integração reduzida 2x2x2, como o C3D20RP/CAX8RP do Abaqus). As funções de
  forma de u são as do próprio mapeamento geométrico (`TPZQuadraticCube::Shape`); as regras são
  `TPZIntCube3D`/`TPZIntQuad`. O espaço H1 de ordem 2 do `TPZCompMesh` não foi usado porque no hexaedro ele
  tem 27 funções (Q2 hierárquico), não as 20 do elemento serendipity do Abaqus/Python.
* Formulação: `[[K_T, -Q], [Qᵀ, S + Δt H]] {Δu, Δp} = -{R_u, R_p}`, Euler implícito + Newton com a
  tangente consistente do Cam-Clay, estado (σ', ε^p, α, ε) por ponto de Gauss, σ'0 e v0 avaliados em cada
  ponto de Gauss (equilíbrio inicial exato) ou no centróide (como as zonas do FLAC3D).
* Sistema linear não simétrico: `TPZSkylNSymMatrix` (LU sem pivotamento; o sistema u-p é quase-definido)
  com os nós renumerados por `TPZCutHillMcKee`.
* Saída VTK legada (`WriteVTK`): `UNSTRUCTURED_GRID` com hexaedros lineares (tipo 12) ou quadráticos de 20 nós
  (`VTK_QUADRATIC_HEXAHEDRON`, tipo 25); campos nodais `deslocamento`, `poropressao`, `excesso_poropressao`,
  `tensao_efetiva`, `tensao_total` (extrapoladas dos pontos de Gauss), `p_efetiva`, `q`; campos por elemento
  `tensao_efetiva_media`, `p_efetiva_media`, `q_medio`, `alpha_medio`, `pc_medio`, `fracao_plastica_passo`,
  `fracao_plastificada`, `poropressao_zona`; `FieldData` com o tempo / fator de carga / δ/H. Os arquivos
  `nome_000.vtk, nome_001.vtk, ...` abrem no ParaView como série temporal.

## Resultados

**Aterro (FLAC3D)**, malha 20x1x10, 10 incrementos de carga não drenados + 25 passos de adensamento até
10⁸ s; saída `aterro_<elem>_000..035.vtk` e `aterro_<elem>_historico.csv`:

* NeoPZ e Python coincidem em todos os instantes (diferenças < 1e-10 m nos recalques e < 1e-8 kPa nas
  poropressões, no limite da precisão dos CSV), com o mesmo equilíbrio global (Σ reações = 4800 kN).
* recalque em x = 0: hex20 0.153 m no fim da etapa não drenada (FLAC3D ≈ 0.14) e 0.275 m em 10⁸ s
  (FLAC3D ≈ 0.19); hex8 0.136 e 0.273 m. A diferença do recalque de adensamento em relação ao FLAC3D já está
  no modelo Python (é dominada pela lei elástica, ver o estudo de sensibilidade em `aterro_itasca.py`), não
  vem do porte.

**Abaqus 1.15.2**, quarto de cilindro 2x2x4 hex20 (321 nós), 150 passos até δ/H = 0.6 em 400 dias; saída
`triaxial_<placa>_<elem>_224_000..015.vtk` e `.csv` (δ/H, p' e q no ponto A, σ_a sob a placa, |p| máx):

* NeoPZ e Python coincidem passo a passo (q, p' no ponto A e σ_a com diferença < 1e-7 kPa e o mesmo
  número de iterações de Newton), para placa lisa e rugosa, hex20 e hex20r.
* placa lisa: q(A) = 149.5 kPa no fim (Abaqus ≈ 148); placa rugosa: q(A) = 141.8 kPa (hex20) e 155.9 kPa
  (hex20r; Abaqus CAX8RP ≈ 153). As tabelas impressas comparam com os pontos digitalizados do Abaqus.

`compara_python.py` refaz os casos com o código Python de `python/` e gera `comparacao_aterro.png` e
`comparacao_triaxial.png`.
