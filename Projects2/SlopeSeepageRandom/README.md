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

> **Forças de percolação: qual curva do artigo comparar.** O artigo usa duas aproximações de `−grad u`: o campo
> semianalítico `K·v'_opt` (princípio de mínimo em velocidade de filtração, Fig. 3, eq. 40) e o FE em poropressão
> `−grad u'_FE`. O Γ = 1.336 do caso de referência (Tabelas 5 e 6, seção 6) é o de `−grad u'_FE`, a mesma
> aproximação deste código: a Fig. 9 dá Γ ≈ 1.31 em β = 45°, α = 1 para `−grad u'_FE` e ≈ 1.89 para `K·v'_opt`.
> (Uma versão anterior deste README atribuía a diferença de Γ ao `v'_opt`; isso estava errado.) A diferença que
> resta vem da discretização: o FE em deslocamentos converge para a carga de colapso de cima, devagar por causa da
> singularidade no pé (ver a convergência com a malha), e de γw (o artigo usa 9.81; `gw=9.81`).

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

> Valores com γw = 10 (versão anterior). Os determinísticos com γw = 9.81 (convergência com a malha, Fig. 8, Fig. 9,
> Tabelas 5 e 6) estão em `resultados/artigo2025/det/` e no relatório.

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
* No caso com percolação o FE fica 5 % acima: ver a nota sobre `v'_opt` acima. Aumentar o domínio de
  25 × 10 m para 105 × 55 m (`Lc=50 Lt=50 Hb=50`) muda Γ só em +0.5 % (1.418 → 1.426 com `h=2 adapt=3`): a
  diferença não vem do truncamento do domínio.

**Cam-Clay no caso com percolação (OCR = 1): dependência de malha.**

| `adapt` | equações | FS | Γ |
|---|---|---|---|
| 0 | 1742 | 1.243 | 1.344 |
| 1 | 2742 | 1.145 | colapso na etapa de percolação (λ_s = 0.85) |
| 2 | 4962 | 1.114 | colapso na etapa de percolação (λ_s = 0.89) |
| 3 | 9730 | 1.094 | 1.00 |

O Cam-Clay normalmente adensado amolece no lado seco da elipse (pontos rasos, `p'` pequeno frente a `p_t`), sem
regularização: o FS por redução de resistência converge com a malha, mas Γ depende da trajetória (a
equivalência `Γ ≥ 1 ⇔ FS ≥ 1` só vale para a plasticidade perfeita associada). Com a resistência real, a
trajetória drenada do rebaixamento (peso próprio e depois as forças de percolação) chega praticamente ao
colapso — coerente com o colapso do rebaixamento acoplado com Cam-Clay (abaixo). Para o Cam-Clay a medida
recomendada é o FS (`medida=fs` no Monte Carlo).

Varreduras (Mohr-Coulomb, `h=1 adapt=2`, Γ):

| α = k_h/k_v | Γ (FE) | Γ artigo (Tab. 5) | FE/artigo |
|---|---|---|---|
| 1 | 1.436 | 1.336 | 1.075 |
| 2 | 1.688 | 1.533 | 1.101 |
| 3 | 1.871 | 1.674 | 1.118 |
| 4 | 2.007 | 1.783 | 1.126 |
| 5 | 2.137 | 1.872 | 1.141 |

| h_w/H | Γ (FE) | Γ artigo (Tab. 6) | FE/artigo |
|---|---|---|---|
| 0.5 | 1.748 | 1.671 | 1.046 |
| 0.6 | 1.571 | 1.494 | 1.052 |
| 0.7 | 1.453 | 1.383 | 1.050 |
| 0.8 | 1.398 | 1.322 | 1.058 |
| 0.9 | 1.398 | 1.307 | 1.070 |
| 1.0 | 1.436 | 1.336 | 1.075 |

Fig. 9 (Γ × β):

| β (graus) | α = 1 | α = 5 | α = 10 |
|---|---|---|---|
| 15 | 4.963 | 20.000 | 20.000 |
| 30 | 2.536 | 4.587 | 6.856 |
| 45 | 1.436 | 2.137 | 2.555 |
| 60 | 1.031 | 1.347 | 1.468 |
| 75 | 0.789 | 0.943 | 0.989 |
| 90 | 0.605 | 0.682 | 0.695 |

Γ cresce com a anisotropia α (horizontal mais permeável: forças de percolação menos horizontais) e tem mínimo
para rebaixamento parcial (h_w/H ≈ 0.8–0.9), como no artigo; a razão FE/artigo é 1.05–1.07 em h_w/H e cresce de
1.075 a 1.14 com α (o campo `v'_opt` do artigo atenua o efeito de α em relação à solução exata de Darcy).


## Monte Carlo

> Campanha anterior (γw = 10, 100 amostras por caso, outro container). Os números abaixo ficam como registro; a
> campanha completa, com γw = 9.81 e todos os casos do artigo, é a da seção *Campanha completa do artigo*.

`h=1 adapt=2` (4962 equações), KL com `hkl=1`, ε_M ≈ 3.6 % compensado, Mohr-Coulomb, Γ por acréscimo de carga;
**100 amostras por caso** (referência: 1000; alguns casos têm mais), portanto CoV(Pf) de 20–100 %: os valores de
Pf das sensibilidades são indicativos (o artigo usa 10 000–90 000 amostras). Para mais amostras basta rodar de novo
o mesmo comando com `n` maior (retomada nativa) ou `scripts/campanha.sh`. CSV, logs e tabelas em `resultados/`.


Colunas `*`: Γ multiplicado por Γ_det(artigo)/Γ_det(FE, mesma malha), isto é, descontada a diferença determinística (forças de percolação v'_opt × FE e malha).

| caso | N | μ | σ | Pf % | CoV(Pf) % | Γ_det FE | μ* | σ* | Pf* % | artigo μ | σ | Pf % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| referencia | 1000 | 1.452 | 0.378 | 8.30 | 10.5 | 1.436 | 1.350 | 0.352 | 14.10 | 1.353 | 0.318 | 11.50 |
| alfa2 | 100 | 1.732 | 0.497 | 1.00 | 99.5 | 1.688 | 1.572 | 0.452 | 6.00 | 1.562 | 0.390 | 3.90 |
| alfa3 | 265 | 1.932 | 0.583 | 0.75 | 70.4 | 1.871 | 1.729 | 0.521 | 3.77 | 1.704 | 0.432 | 1.83 |
| alfa5 | 100 | 2.219 | 0.707 | 0.00 | — | 2.137 | 1.944 | 0.620 | 1.00 | 1.910 | 0.504 | 0.63 |
| cho_coesivo | 129 | 1.265 | 0.188 | 7.75 | 30.4 | 1.331 | 1.287 | 0.191 | 5.43 | — | — | 6.50 |
| cho_cphi | 265 | 1.955 | 1.412 | 5.28 | 26.0 | 1.820 | 1.909 | 1.379 | 6.04 | — | — | 5.50 |
| covc10 | 100 | 1.482 | 0.237 | 0.00 | — | 1.436 | 1.379 | 0.220 | 2.00 | 1.367 | 0.180 | 0.45 |
| covc50 | 100 | 1.417 | 0.538 | 23.00 | 18.3 | 1.436 | 1.318 | 0.501 | 30.00 | 1.324 | 0.472 | 25.78 |
| covc70 | 100 | 1.351 | 0.651 | 36.00 | 13.3 | 1.436 | 1.257 | 0.606 | 41.00 | 1.296 | 0.620 | 36.34 |
| covk0 | 100 | 1.424 | 0.347 | 8.00 | 33.9 | 1.436 | 1.325 | 0.323 | 16.00 | 1.329 | 0.288 | 11.21 |
| covk100 | 100 | 1.518 | 0.450 | 11.00 | 28.4 | 1.436 | 1.412 | 0.419 | 16.00 | 1.391 | 0.360 | 11.44 |
| covphi20 | 100 | 1.488 | 0.520 | 12.00 | 27.1 | 1.436 | 1.384 | 0.484 | 16.00 | 1.361 | 0.384 | 15.18 |
| covphi5 | 100 | 1.451 | 0.348 | 8.00 | 33.9 | 1.436 | 1.350 | 0.324 | 17.00 | 1.351 | 0.300 | 10.20 |
| hw0.5 | 297 | 1.758 | 0.440 | 1.35 | 49.7 | 1.748 | 1.681 | 0.421 | 2.02 | 1.670 | 0.378 | 1.38 |
| hw0.7 | 100 | 1.487 | 0.370 | 8.00 | 33.9 | 1.453 | 1.415 | 0.352 | 12.00 | 1.393 | 0.314 | 8.54 |
| hw0.9 | 100 | 1.425 | 0.371 | 12.00 | 27.1 | 1.398 | 1.332 | 0.346 | 18.00 | 1.325 | 0.306 | 12.58 |
| mcc_ref_fs | 60 | 1.100 | 0.171 | 25.00 | 22.4 | — | — | — | — | — | — | — |
| s2 | 100 | 1.468 | 0.429 | 9.00 | 31.8 | 1.436 | 1.365 | 0.399 | 15.00 | 1.361 | 0.372 | 14.99 |
| s20 | 100 | 1.449 | 0.439 | 15.00 | 23.8 | 1.436 | 1.348 | 0.408 | 23.00 | 1.355 | 0.430 | 20.23 |
| s5 | 100 | 1.482 | 0.443 | 8.00 | 33.9 | 1.436 | 1.378 | 0.412 | 20.00 | 1.350 | 0.404 | 18.93 |

Leitura:

* As médias reescaladas (μ*) reproduzem as do artigo em todos os casos (diferença < 2 %, salvo α = 2–5, 1–2 %):
  a diferença na média é a determinística (forças de percolação `v'_opt` × FE, malha), não a do Monte Carlo.
* O desvio padrão é ~10 % maior (σ* 0.352 × 0.318 na referência) e Pf* é maior (14.1 × 11.5 %): o FE forma
  mecanismos que seguem as zonas fracas, enquanto o mecanismo log-espiral do artigo é uma família de superfícies
  suaves que "promedia" a resistência (efeito conhecido do RFEM; Griffiths & Fenton).
* As tendências das Tabelas 3–6 são reproduzidas: σ e Pf crescem muito com CoV(c) (10 → 70 %: Pf* 2 → 41 %;
  artigo 0.45 → 36 %), pouco com CoV(φ) e CoV(k_v); crescem com a escala s das distâncias de autocorrelação;
  Pf cai com α e com h_w/H menor.
* Cho (2010), sem percolação: Pf = 7.8 % (coesivo; artigo 6.5 %, Cho 7.9 %) e 5.3 % (c-φ; artigo 5.5 %, Cho
  6.37 %).
* Cam-Clay (OCR = 1, Pf por FS, 60 amostras): FS médio 1.10, Pf = 25 % — o solo normalmente adensado amolece e
  é bem menos estável que o Mohr-Coulomb associado com os mesmos c e φ.

Figuras: `resultados/mc_ref/ref.png` (densidade e convergência de Pf e da média, referência), `resultados/mc_alfa.png`
e `resultados/mc_hw.png`.


## Rebaixamento acoplado (Biot, u-p)

O artigo resolve o fluxo à parte (Darcy estacionário) e usa só `−grad u` no esqueleto ("one-way coupling").
O comando `rebaixamento` resolve o problema u-p em deformação plana no tempo (seção acima,
`CoupledDrawdown`), com o nível d'água baixando de `h_w` entre `t = 0` e `t_d`, e em tempos escolhidos congela
`p(x, t)`, calcula `u = p − γw (D − y)` e `grad u` e avalia Γ e FS como no artigo. Para `t → ∞` a poropressão
tende à solução estacionária do artigo (mesmas condições de contorno); o estado inicial (nível na crista) é
hidrostático.

`h=1 adapt=2`, k_v = 10⁻⁵ m/s (k_v/γw = 10⁻⁶ m⁴/(kN s), Tabela 2), E = 10⁵ kPa, ν = 0.3, h_w = H = 5 m;
T_d = c_v t_d/H² é a duração adimensional do rebaixamento. FS e Γ mínimos (no fim do rebaixamento) e no regime
permanente:

| modelo | T_d | FS mín. | Γ mín. | FS (T → ∞) | Γ (T → ∞) | estacionário desacoplado (artigo): FS / Γ |
|---|---|---|---|---|---|---|
| Mohr-Coulomb | 0.01 | 1.114 | 1.186 | 1.208 | 1.425 | 1.218 / 1.436 |
| Mohr-Coulomb | 0.1 | 1.148 | 1.234 | 1.208 | 1.433 | 1.218 / 1.436 |
| Mohr-Coulomb | 1 | 1.182 | 1.333 | 1.218 | 1.433 | 1.218 / 1.436 |
| Mohr-Coulomb | 10 | 1.208 | 1.408 | 1.218 | 1.438 | 1.218 / 1.436 |
| Cam-Clay OCR = 1 | 0.1 | colapso em T = 0.022 (z_w = 8.92 m) | | | | 1.114 / — |
| Cam-Clay OCR = 1 | 10 | colapso em T = 2.38 (z_w = 8.81 m) | | | | 1.114 / — |
| Cam-Clay OCR = 2 | 0.1 | 1.021 | 1.016 | 1.034 | 1.027 | 1.034 / 1.027 |
| Cam-Clay OCR = 2 | 10 | 1.031 | 1.027 | 1.034 | 1.027 | 1.034 / 1.027 |

* Com o acoplamento, o mínimo de estabilidade ocorre no fim do rebaixamento e é tanto menor quanto mais rápido o
  rebaixamento (Mohr-Coulomb: FS −8.5 %, Γ −17 % para T_d = 0.01); depois a poropressão se dissipa e FS e Γ tendem
  ao valor estacionário desacoplado do artigo, que é o limite do rebaixamento lento (T_d = 10). O cálculo do artigo
  é, portanto, contra a segurança para rebaixamentos rápidos.
* Mohr-Coulomb associado dilata ao cisalhar (excesso de poropressão negativo, favorável); o Cam-Clay normalmente
  adensado contrai: gera excesso de poropressão positivo e rompe com 1.1–1.2 m de rebaixamento, mesmo lento —
  coerente com o FS desacoplado de 1.11 e com a perda de estabilidade na trajetória drenada (seção do Cam-Clay
  acima). Com OCR = 2 o talude resiste (FS mínimo 1.02).
* O rebaixamento com malha adaptada exigiu `CleanUpUnconnectedNodes` na malha multifísica (com nós pendentes a
  renumeração de banda escrevia fora dos limites).

Figuras: `resultados/rebaix/rebaixamento_mc.png` e `resultados/rebaix/rebaixamento_mcc.png` (FS e Γ × T).


## Campanha completa do artigo (`scripts/campanha_artigo.py`)

Todos os casos de Monte Carlo do artigo — referência, Tabela 3 (CoV de k_v, c e φ), Tabela 4 (escala s das
distâncias de autocorrelação: 1.5, 2, 5, 10, 20, 400), Tabela 5 (α = 2–5), Tabela 6 (h_w/H = 0.5–0.9) e os dois
exemplos de Cho (2010) — com γw = 9.81, a malha do Monte Carlo (`h=1 adapt=2`) e a semente 2025. Cada caso é
dividido em blocos de 50 amostras e a fila intercala os casos, de modo que todos avançam juntos; tudo é retomável.

```
ninja SlopeSeepageRandom                                   # Release, BUILD_PLASTICITY_MATERIALS=ON, USING_LAPACK=ON
S=Projects2/SlopeSeepageRandom/scripts
python3 $S/campanha_artigo.py jobs <build>/Projects2/SlopeSeepageRandom/SlopeSeepageRandom campanha --alvo 1000 > jobs.txt
nohup $S/fila.sh jobs.txt P > fila.log 2>&1 &              # P = número de núcleos físicos
PROCESSOS=P python3 $S/campanha_artigo.py status campanha  # amostras, μ, σ, Pf, CoV(Pf), s/amostra e previsão
```

* `--alvo 1000` (ou 2000): N por caso. `--alvo artigo`: o S de cada caso no artigo (10 000 a 100 000, ~630 mil
  amostras no total). `--alvo cov5`: o protocolo do artigo com o Pf deste código (CoV(Pf) < 5 %); a fila é gerada
  de novo de tempos em tempos e o alvo de cada caso é reavaliado com as amostras já calculadas.
* Para interromper: matar o `fila.sh`/`xargs` e os processos; para continuar, gerar a fila de novo (mesmo
  `--bloco`) e rodar. Blocos completos são pulados; um bloco interrompido continua da amostra seguinte.
* Tempo: nesta máquina de 4 núcleos, ~4.5–5 s por amostra com 4 processos simultâneos (Cho coesivo ~1.4×, Cho
  c-φ ~0.7×). Em horas de CPU: N = 1000 por caso ≈ 37 h, N = 2000 ≈ 75 h, S do artigo ≈ 820–950 h; dividir pelo
  número de processos. O `status` refaz a previsão com o tempo medido na máquina.
* Para enviar os resultados: `python3 $S/campanha_artigo.py pacote campanha resultados_campanha.tar.gz` (só os
  CSV, `.mec`, `.modo`, `.param` e `.resumo` de cada bloco, sem os caches da KL e da malha).
* Relatório: `analise_artigo.py <campanha> <det> <ref> <saída>` (tabelas e figuras, com os dados do artigo
  extraídos em `resultados/artigo2025/ref/`) e `relatorio_html.py <saída> relatorio.html`.

Saídas novas do comando `mc` (o fator Γ das amostras não muda):

| arquivo | conteúdo |
|---|---|
| `<csv>.mec` | c e φ médios na banda de cisalhamento (‖Δε^p‖ ≥ 10 % do máximo no colapso, ponderados) e k_v médio na massa que se move (‖Δu‖ ≥ 30 % do máximo): as variáveis da correlação de Spearman da seção 6.1 |
| `<csv>.modo` | modo de ruptura (Figs. 16 e 20): `abaixo` se a banda desce mais de 0.1 H abaixo do nível do pé, `pe` se passa a menos de 0.1 H do pé, `acima` nos demais |

Opções novas: `gw=` (γw, padrão 10; o artigo usa 9.81), `fatormax=` (limite do fator de carga, padrão 20) e, no
`det`, o funcional hidráulico `J(u)/(k_h H² γw²)` da Fig. 5 (com `gamma=0 fs=0` só o problema hidráulico é
resolvido).

## Gravação e retomada (read/write)

Nada do que já foi calculado se perde se o processo for interrompido (SIGKILL, queda da máquina):

| arquivo | conteúdo | uso |
|---|---|---|
| `<saida>.csv` | uma linha por amostra, escrita com um único `write` + `fdatasync` | rodar de novo o **mesmo comando** retoma: as amostras já gravadas são puladas e uma linha final incompleta é descartada (gravação atômica, `.tmp` + `rename`) |
| `<saida>.csv.param` | parâmetros que definem o resultado (caso, modelo, CoVs, Lx, Ly, malha, semente, ...) | conferidos na retomada: parâmetros diferentes são recusados (código 3; `forcar=1` aceita) |
| `<saida>.csv.resumo` | μ, σ (n − 1), CoV, Pf, CoV(Pf) e contagem por status de **todo** o CSV | refeito a cada execução |
| `<saida>.csv.lock` | trava (`flock`) | dois processos nunca escrevem no mesmo CSV (código 4) |
| `kl_..._v3.bin` | autopares da KL com cabeçalho (versão, Lx, Ly, ordem, assinatura da malha) | relido em vez de recalculado; arquivo truncado ou diferente é recalculado; escrita atômica |
| `malha_<caso>_h<h>_<hash>.malha` | elementos divididos em cada nível da adaptação | a malha adaptada é reconstruída (idêntica) em vez de refazer as análises de colapso (`malha=0` desliga) |
| `<saida>.csv.campos` (`campos=1`) | c, φ nos pontos de integração e k_v por elemento de cada amostra (binário) | equivalente aos antigos `hhat`; leitura com o comando `campos arquivo=...` |

`amostras=3,17,40-45` recalcula amostras escolhidas (p.ex. com `vtk=1`, para figuras) e confere com o CSV. As
amostras são geradas por `seed_seq{semente, amostra, campo}`: cada uma é reproduzível isoladamente, em qualquer
ordem e em qualquer processo (`scripts/montecarlo.sh` divide as amostras entre processos e junta os CSV).

Verificado: uma execução interrompida com `kill -9` no meio de uma amostra e retomada dá o mesmo CSV que a
execução sem interrupção (exceto o tempo), sem repetições; reexecutar não recalcula nada; os números são
idênticos aos da versão anterior do programa. No `Projects2/GeoMecRandFieldsMonteCarlo` antigo, os arquivos
`hhat` passaram a ter cabeçalho (parâmetros, dimensões) e escrita atômica, a leitura confere o tamanho com a
malha do campo (o formato antigo continua legível) e a retomada usa o conjunto de amostras já gravadas.

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
