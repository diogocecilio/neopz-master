# SlopeMohrCoulomb

Fator de segurança (FS) de um talude em deformação plana com Mohr–Coulomb associativo
perfeitamente plástico, por **aumento de carga (gravidade)** e por **redução de resistência (SRM)**.

* Malha: a mesma do projeto original (`TriGMesh` de `GeoMecDeterm`): domínio 70 × 40 m, talude de 10 m a 45°,
  triângulos P2; base fixa, laterais em rolete.
* Solo: E = 20000 kPa, ν = 0,49, c = 10 kPa, φ = ψ = 30°, γ = 20 kN/m³.
* Modelo constitutivo padrão: projeção analítica em espaço Haigh–Westergaard rotacionado com
  tangente consistente (Lira Cecílio, *J. Eng. Math.* 157:10, 2026) —
  `TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>`. A opção `pv` usa o modelo iterativo antigo
  `TPZPlasticStepPV<TPZYCMohrCoulombPV>` (corrigido) para verificação cruzada.

## Algoritmo

`SlopeAnalysis.h` (único arquivo de driver):

1. `Reset()`: memória plástica virgem em todos os materiais e análise nova.
2. `Newton()`: Newton–Raphson com a tangente consistente, busca linear por retrocesso em ‖R‖,
   critério ‖R‖ ≤ 10⁻⁸ ‖F_ext‖.
3. `Continuation()`: aumenta um parâmetro p a partir de um estado de equilíbrio; se Newton
   converge, aceita o estado (memória plástica) e aumenta o passo; se falha, descarta a tentativa
   e divide o passo por 2. Para quando o passo < 0,2 % de p. O FS é o último p convergido
   (avisa se não houve colapso até p_max).
   * **Aumento de gravidade**: p = λ multiplica γ (resistência intacta).
   * **SRM**: gravidade total aplicada por continuação; depois p = F com c/F e tan φ/F (também ψ).
4. `MarkPlasticZone()` + `Refine()`: marca os elementos com √J₂(εᵖ) ≥ 10 % do máximo no colapso
   do GI e do SRM (mecanismos diferentes) e faz refinamento h com balanceamento 2:1 e renumeração
   de banda (lógica de `DivideElementsAbove`/`Hrefine` do projeto original).
5. `PostPlasticity()` / `CreatePostProcessingMesh()` / `PostProcessVariables()` (as do projeto
   original): VTK com `POrder`, `Atrito` (rad) e `Coesion` efetivos no ponto (com a redução do SRM),
   `StrainPlasticJ2`, `FailureType` (0 elástico, 1 plano principal, 2/3 arestas, −1 ápice) e o
   deslocamento total. Os campos são projetados em L² pelo `TPZPostProcAnalysis` (pequenas
   oscilações em `StrainPlasticJ2`/`FailureType` são da projeção).

## Compilar e executar

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON <neopz>
ninja SlopeMohrCoulomb
./SlopeMohrCoulomb            # modelo do artigo, 5 ciclos de refinamento
./SlopeMohrCoulomb pv         # modelo antigo (verificação)
./SlopeMohrCoulomb check      # teste de Taylor da tangente (Eq. 63-65 do artigo)
./SlopeMohrCoulomb nref=0     # sem refinamento
./SlopeMohrCoulomb nu=0.3     # outro coeficiente de Poisson
```

Saídas: `slope_*_GI_refK.scal_vec.0.vtk`, `slope_*_SRM_refK.scal_vec.0.vtk` (campos no último
estado convergido) e `slope_*_GI_refK.vtk`, `slope_*_SRM_refK.vtk` (malha geométrica).

## Resultados

Teste de Taylor da tangente no ponto material (`check`, Eq. 63–65 do artigo; estados em referencial
girado, inclusive autovalores repetidos):

| Regime | PV2/Voigt (artigo) | PV legado antes | PV legado corrigido |
|---|---|---|---|
| elástico | exato | ordem 1,00 | exato |
| plano principal | 2,00 | 1,00 | 2,00 |
| aresta σ₂ = σ₃ | 2,00 | 1,00 | 2,00 |
| aresta σ₁ = σ₂ | 2,00 | 1,00 | 2,00 |
| ápice | exato (Dep = 0) | exato | exato |

Tangentes simétricas (|D−Dᵀ|/|D| < 2·10⁻¹⁴), o que legitima a skyline LDLᵀ. Na análise global o
Newton converge quadraticamente (p.ex. ‖R‖/‖F‖: 4·10⁻⁵ → 4·10⁻⁸ → 1·10⁻¹³).

Fator de segurança (ν = 0,49, tolerância de 0,2 % no parâmetro; `nref=5`):

| ciclo | equações | FS gravidade (artigo) | FS SRM (artigo) | FS gravidade (PV) | FS SRM (PV) |
|---|---|---|---|---|---|
| 0 | 870 | 3,043 | 1,396 | 3,043 | 1,401 |
| 1 | 918 | 2,516 | 1,312 | 2,516 | 1,312 |
| 2 | 1140 | 2,094 | 1,258 | 2,094 | 1,258 |
| 3 | 1814 | 1,918 | 1,229 | 1,918 | 1,229 |
| 4 | 3408 | 1,840 | 1,215 | 1,840 | 1,212 |
| 5 | 7684 / 7142 | 1,793 | 1,203 | 1,797 | 1,203 |
| Bishop simplificado | — | 1,857 | 1,209 | | |

* Os dois modelos constitutivos, independentes, dão o mesmo FS (diferenças dentro da tolerância).
* O FS converge com o refinamento da zona plástica para os valores de equilíbrio-limite (Bishop,
  busca de círculos, script independente). Malha grossa superestima muito o FS por gravidade
  (mecanismo raso, junto à face).
* O fator de gravidade não é o FS de redução de resistência: para Mohr–Coulomb, multiplicar γ por λ
  equivale a dividir só c por λ, logo FS_gravidade ≥ FS_SRM quando φ > 0.
* Sensibilidade ao coeficiente de Poisson (`nu=0.3`, modelo do artigo, ciclo 5): FS gravidade 1,789
  e FS SRM 1,207 (contra 1,793 e 1,203 com ν = 0,49): sem travamento volumétrico relevante.
* Incerteza do FS pela continuação: < 0,4 % (último passo convergido e falha a ≤ 2·0,2 %).

## Bugs e inconsistências corrigidos na biblioteca

Convenção adotada (a do artigo e de `TPZElasticResponse::De`): deformações em Voigt com
cisalhamento de engenharia (γ = 2ε) nos componentes XY, XZ, YZ do `TPZTensor`; tensões com
componentes tensoriais. O commit `e16d8bc` introduziu essa convenção no material e em
`TPZElasticResponse`, mas não em `TPZPlasticStepPV` nem no pós-processamento.

| Arquivo | Problema | Efeito | Correção |
|---|---|---|---|
| `TPZPlasticStepPV.cpp` | Autovetores de ε_tr (com γ no slot XY) e tangente derivada em relação a ε tensorial | tangente errada **até no regime elástico** (ordem de Taylor 1, D(XY,XY)=2G), Newton linear | autossistema único do σ de teste (εᵢ−εⱼ = (σᵢ−σⱼ)/2G) e colunas de cisalhamento ×½ (d/dγ) |
| `TPZPlasticStepPV.cpp` | Termo de rotação da tangente com os autovalores **de teste** de σ (deveriam ser os projetados, Souza Neto C.14) | rigidez de rotação elástica em todo ponto plástico | autovalores projetados no `TangentOperator` |
| `TPZPlasticStepPV.cpp` | `tempMat` acumulado entre pares (i,j) e `Tangent +=` sem zerar | tangente 3D errada; lixo de memória em `Dep(6,6)` não inicializado | matriz por par; `Tangent.Redim(6,6)` |
| `TPZPlasticStepPV.h/.cpp` | `fReductionFactor()` = 0; SRM só aplicado se `fmatprop[0] > 1e-3` | c/0 → material elástico; sem propriedades por ponto a redução era **ignorada** | padrão 1; `LocalCriterion()`: cópia do critério com as propriedades do ponto (ignora vetor nulo de placeholder) e a redução — nunca compõe |
| `TPZPlasticStepPV.cpp` | `ApplyStrainComputeDep` chamava `ProjectSigmaDep` (`DebugStop` no Mohr–Coulomb); `Phi()` sem parâmetros locais | `TaylorCheck` inutilizável para MC; `Yield` pós-processado com resistência não reduzida | especialização para MC via `ApplyStrainComputeSigma` (demais critérios inalterados); `Phi` usa `LocalCriterion()` |
| `TPZYCMohrCoulombPV.h/.cpp` | `GetLocalMatState` com `DebugStop`; `ChangeLocalMatParameters` recomputava de `fmatprop` (sem propriedades por ponto não reduzia nada) | SRM ignorado sem `fmatprop` | redução c/F, atan(tanφ/F), atan(tanψ/F) dos parâmetros correntes, aplicada a uma cópia; `SetLocalMatState` mantém ψ = φ (associativo, como antes) |
| `TPZYCMohrCoulombPV.cpp` | Validade das arestas com `IsZero` absoluto (1e-12) | arredondamento manda retorno de aresta válido ao ápice (tração) | tolerância relativa ao nível de tensão |
| `TPZYCMohrCoulombPV.cpp` | ψ = 0 no ápice: 0/0; `k_proj` não definido no ramo elástico; `SigmaElastPV` atribuía o vetor inteiro | NaN; endurecimento zerado | guardas; `k_proj = k_prev`; índices |
| `TPZYCMohrCoulombPV2.cpp` | Ramo elástico `alphan = alphan1` (sobrescreve a **entrada** com a saída não inicializada) | memória com lixo (UB) | `alphan1 = alphan` |
| `TPZYCMohrCoulombPV2.cpp` | Ápice com Jacobiano −1e-11 (artigo: Dep = 0); teste elástico `Φ < 0` (artigo Eq. 56: `Φ ≤ 0`); todas as regiões com `m_type = 1` | tangente negativa espúria; regiões indistinguíveis | Jacobiano nulo; `≤ 0`; 1 plano, 2 aresta direita, 3 aresta esquerda, −1 ápice |
| `TPZYCMohrCoulombPV2.h/.cpp` | `YieldFunction`, `SetLocalMatState`, `ChangeLocalMatParameters`, `Print` vazios; funções não-void sem `return`; ψ ignorado sem aviso | SRM impossível; UB | implementados (Eq. 44, 45, 49; c/F, atan(tanφ/F)); stubs com `DebugStop`; `DebugStop` se ψ ≠ φ (método associativo) |
| `TPZPlasticStepVoigt.cpp/.h` | `STATE hvarnew;` não inicializado gravado em `fN.m_hardening` | UB / NaN na memória | inicializado com o estado anterior |
| `TPZPlasticStepVoigt.cpp` | `SetElasticResponse` não repassava o ER ao critério (G, K da projeção) | projeção e tangente com módulos diferentes do preditor (ordem 1) | repassa (também no construtor) |
| `TPZPlasticStepVoigt.cpp/.h` | Sem fator de redução; cópia descartava `fN`; `ClassId()==0`, `Read/Write` vazios; `Phi()`=0; `ApplyStrainComputeDep` nulo | SRM impossível; serialização e `Yield` inúteis | `SetStrengthReductionFactor`, `LocalCriterion()`, cópia completa, serialização, `Phi` real |
| `TPZPlasticStepVoigt.cpp` | Ramo elástico reconstruía σ pela decomposição espectral | erro ~1e-8 com autovalores repetidos (atrapalha o critério de Newton) | σ = σ_trial; no ramo plástico soma só a correção plástica |
| `TPZYCTrescaVoigt.cpp/.h`, `TPZYCVonMisesVoigt.h` | mesmo `alphan = alphan1` do PV2 (Tresca); `YieldFunction` = `DebugStop` | endurecimento corrompido; `Phi()`/`Yield` abortaria | `alphan1 = alphan`; funções de escoamento implementadas |
| `TPZElasticResponse.cpp` | `operator=` não copiava σ* | cópia incompleta | corrigido |
| `TPZPlasticState.h` | Construtor de cópia: `fmatpropinit(source.fmatprop)` | propriedades "iniciais" viram as reduzidas → reduções compostas | `fmatpropinit(source.fmatpropinit)` |
| `TPZMatElastoPlastic_impl.h` | `m_m_type` não gravado na memória; `m_u` guardava só o último incremento | `FailureType` sempre 0; `DisplacementDoF` errado | grava `m_m_type`; `m_u += Δu` |
| `TPZMatElastoPlastic_impl.h` | Elasticidade não linear gravava `m_ER` na memória a cada avaliação (fora de `fUpdateMem`) | escrita fora do commit (sem efeito observável hoje) | ER local; gravado só no commit |
| `TPZMatElastoPlastic_impl.h` | Autovalores e J₂ de deformação calculados com γ no slot tensorial | J₂(εᵖ) com parcela de cisalhamento ×4 | conversão γ→ε antes dos invariantes |
| `TPZMatElastoPlastic_impl.h/.h` | Cópia sem `m_force0`/`fExactSolution`; ponteiro não inicializado; `EEXACT` lia `sol[1]`; `NSolutionVariables` 6 para 3 valores | UB, pós-processamento errado | corrigidos |
| `TPZMatElastoPlastic_impl.h/.h`, `TPZYCMohrCoulombPV2.h` | `POrder`, `Coesion`, `Atrito` pedidos pelo pós-processamento original não existiam | variável desconhecida | `EPOrder`, `ECohesion`, `EFriction` (c e φ do critério local, já reduzidos) |
| `TPZMatElastoPlastic2D_impl.h` | `auto bc_with_memory = ...` (cópia profunda de toda a memória do contorno a cada ponto de integração); BC tipo 4 escrevia em `Val2` compartilhado | custo O(n²) por montagem; corrida entre threads | referência; tração local |
| `TPZMatWithMem.h` | `operator=` sem `return` e chamada inválida em `shared_ptr` | UB / não compila se usado | corrigido |
| `pzelastoplasticanalysis.cpp` | Buscas lineares dicotômica (padrão) e áurea clonavam **a malha inteira com a memória** a cada iteração e vazavam; a dicotômica nunca testava o passo completo | vazamento de memória; convergência linear | avaliação in loco; na dicotômica o passo completo é aceito se reduzir ‖R‖ |
| `pzelastoplasticanalysis.cpp` | `ELineSearch::None` deixava `fSolution` vazio; sem saída para NaN (agora mantém o último iterado finito); `ManageIterativeProcess` aceitava passo não convergido; `SetAllCreateFunctionsWithMem` instalava 8 ponteiros nulos | falhas/estado plástico corrompido | corrigidos |
| `pznonlinanalysis.cpp` | Regressão: `sol += prevsol` em cópia local descartada | Newton não-linear sem busca linear nunca converge | `fSolution += prevsol` |

## Artigo × código (Lira Cecílio, 2026)

Verificado numericamente (transcrição literal em numpy, referência independente de Souza Neto
Box 8.5 e projeção de ponto mais próximo exata por conjunto ativo; 24 000 estados aleatórios):

* Tabelas 4–6, Eq. 57–58 e o código C++ coincidem (tensão 1e-14, Jacobiano 4e-16); seleção de
  região (Eq. 59–61) idêntica; tangente (Eq. 29–37 e Listing 1) idêntica ao código e às diferenças
  finitas (1,7e-8); teste de Taylor com ordem 2,00 no plano principal e nas duas arestas.
* Divergências código × artigo (corrigidas no código): Jacobiano do ápice −1e-11 (artigo: Dep = 0);
  teste elástico `Φ₁ < 0` (artigo: `Φ₁ ≤ 0`); ψ ignorado sem aviso (o método é associativo —
  agora `DebugStop` se ψ ≠ φ). Sem efeito: tolerância ε_κ absoluta 1e-12 no C++ e 1e-15 no
  Listing 1 (nas arestas κ = 0 por ambos os ramos).
* Exemplo da sapata (`Projects2/PlasticityTestsMohrCoulomb`, removido do branch): o driver aplicava pressão
  (Neumann, tipo 1) e não o recalque prescrito da Seção 5.1; não calculava P = R/B; a regra de
  integração do P2 no NeoPZ é a de 6 pontos (ordem 2p), o artigo usa 3 pontos; o penalty
  `fBigNumber` ≈ 6,7·10¹⁶ com deslocamento prescrito e tolerâncias absolutas de 1e-6 em
  `NewtonRaphson` pode impedir a convergência; aceitava incrementos não convergidos.
* Erratas no artigo (não afetam o código):
  * Eq. (54) está multiplicada por 432 (falta o fator 1/432); minimizador inalterado.
  * Tabelas 4–6 trazem "sin φ²" onde deve ser sin²φ.
  * Eq. (28) deveria usar (v_jj^ε)ᵀ (mapa de deformação) e é apenas a parte coaxial; a Eq. (29)
    implementada (A + R_V) está correta.

## Driver original (`Projects2/GeoMecDeterm`, arquivos enviados) — por que foi reescrito

* SRM inalcançável: `LoadingRamp`/`ApplyGravityLoad` eram `DebugStop`; sem fmatprop a redução
  era ignorada pela biblioteca; `ShearRed` gravava `fmatprop = c₀/FS` e a biblioteca dividia de
  novo (FS²).
* Dentro de cada bissecção as tentativas não eram independentes: `InitializeMemory` não zera εᵖ/σ e
  cada tentativa partia do último estado aceito.
* Critérios absolutos (‖R‖ < 1e-6) dependentes de unidades e da malha; sem detecção de NaN; o
  fallback reiniciava do iterado divergido.
* `FS`/`FSOLD` não inicializados em `Solve`/`SolveDeterministic` (`FSOLD` vira `FS0` de `ShearRed*`;
  se 0, `FS = 1/((1/FSmin + 1/FSmax)/2)` divide por zero); arc-length aceitava passos não
  convergidos e o teste |Δλ| < tol parava antes do pico.
* Construtor de cópia raso + destrutor que apaga as malhas (double free); as variáveis de
  pós-processamento `Atrito`, `Coesion` e `POrder` não existiam no material (agora implementadas);
  deslocamentos sempre nulos após `AcceptSolution` (agora o pós-processamento usa a solução
  acumulada). `PostPlasticity`/`CreatePostProcessingMesh`/`PostProcessVariables` foram mantidas.
* Skyline simétrica + LDLᵀ com a tangente não simétrica do caminho PV (descartava o triângulo
  inferior); sem renumeração de banda, a skyline explodia após o refinamento.

## Pendências conhecidas (documentadas, não alteradas)

* `TPZPlasticStepPV::ApplyStressComputeStrain` e `ApplyLoad` (não usados) continuam com lógica
  incorreta; `TPZPlasticStepPV::Write/Read` não serializam o fator de redução (agora padrão 1).
* `TPZPlasticStepPV` com CamClay/DruckerPrager/Sandler mantém o caminho antigo da tangente
  (`ProjectSigmaDep`), que não segue a convenção de cisalhamento de engenharia.
* `e16d8bc` fez `TPZElasticResponse::ComputeStress/ComputeStrain` ignorarem ε*/σ*. Restaurar exige
  corrigir também `TPZPorousElasticResponse::ComputeStress` (usa μ·dev(ε) nas normais, deveria ser
  2μ), senão a elasticidade porosa regride; não alterado.
