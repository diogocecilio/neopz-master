# Talude com forças de percolação e variabilidade espacial (NeoPZ)

Reprodução, por elementos finitos elastoplásticos no NeoPZ, de

> M. Vargas Ceron, D. L. Cecílio, R. V. Linn, S. Maghous, *Stability Analysis of Slope Subjected to Seepage
> Forces Considering Spatial Variability of Soil Properties*, Int. J. Numer. Anal. Methods Geomech. 49 (2025)
> 2459–2491,

com **Mohr-Coulomb** e **Cam-Clay modificado**, campos aleatórios de Karhunen-Loève (rotinas do
`Projects2/GeoMecRandFieldsMonteCarlo`, `ElasticityRandomField` e `eigensolver`, reescritas com as classes
nativas) e, além do artigo, o **rebaixamento acoplado** (Biot, u-p) com o material
`TPZMatPoroElastoPlastic3DMem` em deformação plana.

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON -DUSING_LAPACK=ON <neopz>     # LAPACK: autoproblema KL (dsbgv)
ninja SlopeSeepageRandom
./SlopeSeepageRandom det caso=percolacao modelo=mc h=1 adapt=3            # Γ e FS determinísticos
./SlopeSeepageRandom mc  caso=percolacao modelo=mc h=1 adapt=2 n=1000      # Monte Carlo (CSV retomável)
./SlopeSeepageRandom rebaixamento modelo=mc h=1 adapt=2 Td=0.1 gamma=1   # u-p acoplado
```

Scripts (`scripts/`): `deterministico.sh` (casos do artigo, Fig. 9 e colunas determinísticas das Tabelas 5 e
6), `montecarlo.sh` (divide as amostras entre processos, retoma e junta os CSV), `rebaixamento.sh` (varredura de
T_d), `fila.sh` (fila de comandos em paralelo) e `analisa_mc.py` (μ, σ, CoV, Pf, CoV(Pf), densidade e
convergência, comparados com as Tabelas 3–6).

## O que é calculado

O artigo avalia a estabilidade pela **análise limite cinemática** (mecanismo log-espiral, limite superior):

* **Γ** (fator de estabilidade, eq. 54): multiplicador das cargas `λ (γ' e_y − grad u)` no colapso;
* **FS** (eq. 58): divisor de `c` e `tan φ` que leva ao colapso com as cargas reais.

`Γ ≥ 1 ⇔ FS ≥ 1`, de modo que `Pf = P(Γ < 1) = P(FS < 1)`; Γ e FS só coincidem para φ = 0. Aqui as duas medidas
são obtidas por **elastoplasticidade incremental** (o colapso é a não convergência do Newton, critério de
Griffiths & Lane, com bissecção do passo até a tolerância `reltol`): Γ por acréscimo das cargas a partir do
estado nulo e FS por redução de resistência (`SetStrengthReductionFactor`, nativo). Para o Mohr-Coulomb
associado a carga de colapso é única (teoremas da análise limite) e o FE converge para o valor exato com o
refinamento, que é menor ou igual ao limite superior do artigo.

**Forças de percolação.** `u = p − γw y'` (excesso de poropressão em relação à hidrostática com o nível na
crista, `y'` = profundidade abaixo da crista) é a solução de Darcy estacionário com `K = k_v e_y⊗e_y +
k_h (1 − e_y⊗e_y)`, `u = 0` na crista, `u = −γw min(y', h_w)` na face e `u = −γw h_w` no pé (eq. 21); base e
laterais impermeáveis. No esqueleto entram `γ'` e `−grad u` como forças de corpo (eq. 5).

> **Diferença importante em relação ao Γ = 1.336 do artigo.** No artigo, `−grad u` usado na análise limite vem
> do campo de velocidades *semianalítico* `v'_opt` (princípio de mínimo em velocidade de filtração, Fig. 3,
> eq. 40, domínio com `R_e`, `L_m = 10 H`); o próprio artigo mostra (eq. 28, Fig. 5, Fig. 9) que esse campo e o
> FE em poropressão dão estimativas *superior e inferior* das forças de percolação. Aqui `u` é o FE em
> poropressão (H1, ordem 2) — a estimativa inferior —, o que leva a Γ maior que o do artigo no caso com
> percolação; nos casos secos (Cho) não há essa diferença e o acordo é de 0.1–2.5 % (abaixo).

## Estrutura (classes nativas do NeoPZ)

| arquivo | conteúdo |
|---|---|
| `SlopeGeometry` | `TPZGeoMesh` estruturado (quadriláteros ou triângulos) do talude, ids de contorno, `Refine`/`Grow` (refinamento h com `TPZGeoEl::Divide` e nós pendentes) |
| `KLRandomField` | KL de Galerkin: `TPZMatKLKernel` + `pzdoublestrmatriz` (matrizes C e B), `TPZSBMatrix` + `TPZLapackEigenSolver::SolveGeneralisedEigenProblem` (dsbgv), modos `Φ_k = √λ_k φ_k`, erro de truncamento `ε_M = 1 − Σλ/|Ω|`, compensação pontual da variância `H/√v(x)`, avaliação por `LoadSolution` + `ComputeSolution`, cache binário (`TPZBFileStream`), lognormal |
| `SeepageProblem` | Darcy anisotrópico (`TPZDarcyFlow` com K por elemento: campo aleatório de k_v), condições de Dirichlet por `SetForcingFunctionBC`, `TPZLinearAnalysis` + skyline LDLᵀ |
| `SlopeStability` | `TPZMatElastoPlastic2D<T, TPZElastoPlasticMem>` + forças de percolação na memória (`fdPorePressure`); `T` = `TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>` (tangente consistente exata) ou `TPZModifiedCamClay`; propriedades por ponto em `TPZPlasticState::fmatprop`; Newton com critério relativo, aceitação por `SetUpdateMem` + `AssembleResidual`; Γ e FS com cortes de passo; `TPZPostProcAnalysis` para VTK |
| `CoupledDrawdown` | u-p acoplado: `TPZMultiphysicsCompMesh` (u H1 vetorial com memória + p H1 linear), `TPZMatPoroElastoPlastic3DMem<T>(id, 2)`, Euler implícito, rebaixamento `z_w(t)` por `SetForcingFunctionBC` |
| `main.cpp` | comandos `det`, `mc` e `rebaixamento`, malha adaptativa, Monte Carlo |

**Malha adaptativa** (`adapt=N`): resolve o problema médio, marca os elementos com `||Δε^p||` do último passo
aceito (o mecanismo de colapso) acima de `frac · max`, acrescenta `camadas` vizinhos e divide (`Divide`). O pé do
talude é um canto reentrante (225°): tensões e gradiente de `u` são singulares ali (`∇u ~ r^−0.2`), e a
convergência de Γ com h é lenta; o refinamento guiado pelo mecanismo é muito mais eficiente que o uniforme.

**Cam-Clay modificado.** Mesma resistência de estado crítico: `M = √3 sin φ` (deformação plana, ou
`6 sin φ/(3 − sin φ)` com `mapeamento=triaxial`), `p_t = c cot φ` (invariante no SRF), tensão inicial de uma
análise elástica geostática na mesma malha, `p_c0` = elipse normalmente adensada por `σ'0` vezes `OCR`, e
elasticidade linear. Com ele Γ e FS dependem da trajetória (endurecimento/amolecimento).

## Resultados determinísticos

`scripts/deterministico.sh` (Mohr-Coulomb, p = 2, malha adaptada `adapt=3`):

| caso | Γ (FE) | Γ artigo | FS (FE) | FS referência |
|---|---|---|---|---|
| Cho (2010) coesivo 2:1, c_u = 23 kPa, H = 5 m | 1.321 | 1.354 | 1.324 | 1.356 (Cho, LE) |
| Cho (2010) c-φ 1:1, c = 10, φ = 30°, **H = 10 m** | 1.775 | 1.777 | 1.207 | 1.204 (Cho), 1.203 (artigo) |
| referência com percolação (Tabela 2), h_w = H = 5 m | 1.402 | 1.336 | 1.202 | — |
| Cho c-φ, Cam-Clay (OCR = 1) | 1.100 | — | 1.132 | — |

* O exemplo c-φ de Cho (2010) tem H = 10 m (o texto da seção 5.3.2 diz 5 m, mas Γ = 1.777 e FS = 1.203/1.204 só
  são reproduzidos com H = 10 m; com H = 5 m o próprio Bishop simplificado dá FS ≈ 1.61).
* No caso coesivo o FE fica 2.4 % abaixo do limite superior log-espiral (esperado).
* No caso com percolação o FE fica 5 % acima: ver a nota sobre `v'_opt` acima.

_(resultados em execução; as tabelas são geradas por `scripts/tabelas.py <diretório de resultados>`)_

## Monte Carlo

_(resultados em execução; as tabelas são geradas por `scripts/tabelas.py <diretório de resultados>`)_

## Rebaixamento acoplado (Biot, u-p)

O artigo resolve o fluxo à parte (Darcy estacionário) e usa só `−grad u` no esqueleto ("one-way coupling").
O comando `rebaixamento` resolve o problema u-p em deformação plana no tempo (seção acima,
`CoupledDrawdown`), com o nível d'água baixando de `h_w` entre `t = 0` e `t_d`, e em tempos escolhidos congela
`p(x, t)`, calcula `u = p − γw (D − y)` e `grad u` e avalia Γ e FS como no artigo. Para `t → ∞` a poropressão
tende à solução estacionária do artigo (mesmas condições de contorno); o estado inicial (nível na crista) é
hidrostático.

_(resultados em execução; as tabelas são geradas por `scripts/tabelas.py <diretório de resultados>`)_

## Correções feitas no NeoPZ e nas rotinas antigas

Biblioteca (usadas por este projeto e pelas rotinas de campos aleatórios):

* `Common/TPZLapack.h`: `USING_LAPACK` não compilava com o LAPACK do sistema (≥ 3.9.1: comprimentos ocultos de
  string; `cblas.h` fora do MKL) — macro `PZ_LAPACK(nome)` em `pzfmatrix`, `pzbndmat`, `pzsbndmat` e
  `TPZLapackEigenSolver`.
* `TPZLapackEigenSolver` (dsbgv): teste de dimensão com `||` trocado por `&&` e `ldbb` lido de A.
* `TPZEigenAnalysis::Assemble`: a matriz C era montada duas vezes (e a primeira vazava).
* `TPZMatKLKernel`: dados de quadratura calculados uma vez por par de elementos (antes `ComputeRequiredData` do
  segundo elemento em cada par de pontos), regra aumentada no bloco diagonal (quina do kernel exponencial),
  regras clonadas sem vazamento; `pzdoublestrmatriz` respeita a opção Galerkin/Nyström. A área de Ω para
  `ε_M` é a integral de 1 (a soma de B não vale para a base hierárquica).
* `TPZPlasticStepVoigt` + `TPZYCMohrCoulombPV2` (o modelo do Monte Carlo antigo): **c e φ por ponto nunca eram
  lidos** de `fmatprop` (todas as amostras resolviam o problema médio) e não havia redução de resistência;
  `hardening` não inicializado no ramo elástico; funções sem `return`; tangente consistente reescrita sem
  matrizes temporárias (mesma fórmula, verificada por diferenças finitas, ~4× mais rápida).
* `TPZPlasticStepPV`: fator de redução padrão 0 (c/0 assim que `fmatprop` era usado) → 1. A tangente desse
  modelo não é consistente com a convenção de distorções de engenharia (o Newton estagna); por isso o projeto
  usa o `TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>`.
* `TPZPlasticState`: o construtor de cópia copiava `fmatprop` em `fmatpropinit` (reduções compostas, FS²);
  `TPZElastoPlasticMem`: a cópia perdia poropressão e gradiente.
* `TPZMatElastoPlastic2D`: `Contribute` sem matrizes temporárias; `ContributeBC` copiava toda a memória do
  contorno em cada ponto de integração.
* `TPZModifiedCamClay`: parâmetros por ponto (c, φ → M, p_t; p_c0; σ0) via `fmatprop` e fator de redução.
* `TPZMatPoroElastoPlastic3DMem`: deformação plana (`dim = 2`), permeabilidade ortótropa, valores de contorno
  por `SetForcingFunctionBC` (nível d'água variável) e instância para o Mohr-Coulomb.

Rotinas antigas do `Projects2/GeoMecRandFieldsMonteCarlo` (substituídas por este projeto): além do
`fmatprop` ignorado, o índice inicial de `FindElement` não era inicializado (`elidsrc`), o pós-processamento
final usava o FS do nível anterior (`FSOLD`), `Val2()[2]` era escrito num vetor de tamanho 2, amostras com
FS > 10 eram descartadas do CSV (Pf com denominador errado), a bissecção devolvia `lo` sem testá-lo, e a
variância truncada da KL (M = 100: 74 % da variância) não era compensada. Aqui: modos até `ε_M` dado,
compensação `H/√v(x)` e todas as amostras gravadas com o status.
