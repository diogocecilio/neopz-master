# SlopeSeepageForces

Estabilidade de taludes sob as forças de percolação de um rebaixamento rápido — reprodução da Seção 4.3 de
Ceron, Cecílio, Linn & Maghous, *Stability analysis of slope subjected to seepage forces considering spatial
variability of soil properties*, IJNAMG 49(11), 2025 (doi:10.1002/nag.3993): Fig. 5 (funcionais hidráulicos),
Fig. 8 (altura crítica × h_w/H) e Fig. 9 (fator de estabilidade × β para α = 1, 5, 10).

Este diretório contém a parte C++ (NeoPZ), construída sobre `Projects/SlopeDrawdown` e
`Projects/SlopeMohrCoulomb`; as implementações de referência em Python estão em `scripts/` (dados digitalizados
em `data/`, resultados em `results/`). O campo analítico K⁻¹·v'_opt (`scripts/analytical_seepage.py`) está
portado (`AnalyticalSeepage.h`, item 5 do *Modelo*); a análise limite (`scripts/limit_analysis.py`) também
(`LimitAnalysis.h`, item 6 do *Modelo*; comandos `la` e `labatch`). As Figs. 5, 8 e 9 são reproduzidas em C++ pelos
comandos `fig5 out=`, `fig8` e `fig9` (CSVs retomáveis em `results/cpp/`, figuras `results/fig5.png`, `fig8.png`,
`fig9.png`; seção *(f)* dos resultados), e verificadas de forma independente pelo MEF elastoplástico (aumento de
gravidade, comando `fembatch`; seção *(g)*).

## Convenções

* Coordenadas NeoPZ (y para cima): aresta do topo O = (0, 0), topo y = 0 para x ≤ 0, face de O ao pé
  T = (H/tan β, −H), terreno do pé y = −H; nível d'água após o rebaixamento no ponto W = (h_w/tan β, −h_w) da face.
  Coordenadas do artigo (y para baixo): x_artigo = x, y_artigo = −y.
* Ids de contorno (os de `TriGMesh`): −1 base, −2 lateral direita, −3 terreno do pé, −4 topo, −5 lateral
  esquerda, −6 face; solo 1.
* u = poropressão em excesso (Eq. 2), p = u − γw y = poropressão total. Força de percolação f = −∇u (kN/m³).

## Modelo

1. **Malhas** (`SlopeGeometry.h`, `DelaunayMesher.h`): triângulos lineares de um polígono com nós em O, W e T,
   por refinamento de Delaunay (Ruppert, ângulo mínimo 28°) com função de tamanho graduada para O, W, T e a face
   (h0, hs, grade, hmax em unidades de H) e refinamento uniforme opcional. Dois usos: domínio hidráulico grande
   (artigo, Fig. 4: 50 H à esquerda de O, 10 H à direita de T, 30 H abaixo de T) e domínio de estabilidade
   (`sa` H + H/tan β de cada lado de O/T e abaixo de T). `SlopeMohrCoulombGMesh` é a `TriGMesh` de
   `SlopeMohrCoulomb` transladada (regressão).
2. **Percolação permanente anisotrópica** (`AnisotropicDarcy.h`, `SeepageFE.h`), Eqs. 20–21: −div(K ∇u) = 0,
   K = diag(k_h, k_v), α = k_h/k_v; H1 de ordem 2 (o `TPZDarcyFlow` do NeoPZ é isotrópico: material próprio
   `TPZAnisotropicDarcy`, Dirichlet por penalidade relativa a K). Dirichlet na superfície (−3, −4, −6):
   u = γw·max(y, −h_w) (u = 0 no topo, u = γw y na face acima d'água, −γw h_w abaixo e no terreno do pé).
   Laterais e base (`hbc`, predefinições de `scripts/fe_seepage.py`): **padrão `zero_lb`** — u = 0 (campo
   distante ainda hidrostático no nível do topo) na lateral esquerda e na base, fluxo nulo na lateral direita.
   É a condição que reproduz as curvas tracejadas da Fig. 5 (0,01–0,22 % em 24 pontos); com o domínio todo
   impermeável J fica 3–9 % abaixo (ver *Resultados*). Saídas: J(u) = ½∫∇u·K∇u (regra dos pontos médios das
   arestas, exata para P2), J/(k_h H² γw²), mínimo de p (deve ser ≥ 0).
3. **Campo de forças de percolação** (`SeepageForceField.h`): `PoreField` (generalização do `PoreField` de
   `SlopeDrawdown`) guarda por triângulo o interpolante quadrático da solução (exato até P2, gradiente P2 exato)
   e localiza pontos por uma quadtree das caixas dos triângulos (~0,2 µs por ponto, sem estado mutável: seguro
   entre threads); fora da malha hidráulica f = 0. Interface genérica
   `slope::ForceField = std::function<void(const TPZVec<REAL>&, REAL f[2])>`, adaptador
   `FromPaperCoordinates` para campos escritos em coordenadas do artigo e `NoSeepage()` (h_w = 0 / seco).
4. **Estabilidade por elementos finitos** (`FEMStability.h`): aumento de gravidade de
   `SlopeMohrCoulomb/SlopeAnalysis.h` (**sem alteração**) na malha de estabilidade com força de volume
   b = λ(γ_ref g + f(x)), g = (0, −1), contorno efetivo livre de tensões. Forma do artigo (`form=u`):
   γ_ref = γ' = γ − γw, f = −∇u; Γ_FEM = λ_crit e H_crit = λ_crit H. Forma de `SlopeDrawdown`
   (`form=p`): γ_ref = γ_sat, f = −∇p, idêntica (b = λ(γ_sat g − ∇p) = λ(γ' g − ∇u)); `form=p+` usa
   p⁺ = max(p, 0) como `SlopeDrawdown`. Como λ chega à função de forçamento: em
   `TPZMatElastoPlastic2D::Contribute` a força local parte de `m_force` (que `SlopeAnalysis` faz λ(0, −γ_ref, 0))
   e a função de forçamento a sobrescreve; ela recupera λ = −m_force[1]/γ_ref (como `SetSeepageForce` de
   `SlopeDrawdown`). Ciclos de refinamento da zona plástica como em `FactorOfSafety` de `SlopeDrawdown`
   (`nref`, `mark`; `srm=1` calcula também o SRM e marca as duas zonas, como lá). `ExtrapolateCycles` extrapola os
   λ dos ciclos para h → 0 (ordem 1, 2 λ_k − λ_{k−1}; Richardson com a ordem observada; reta em 1/√neq; ver *(g)*).
5. **Campo analítico K⁻¹·v'_opt** (`AnalyticalSeepage.h`, porte de `scripts/analytical_seepage.py`; dedução e
   correções das Eqs. 31 e 39 impressas em `scripts/analytical_seepage_derivation.md`): velocidade ótima da classe
   da Eq. 29 em torno de O — zona 1 (r < R_w, Eq. 31 corrigida: prefator 1/(D − C), expoente e coeficiente √(C/D)
   no termo e_θ), zona 2 (R_w ≤ r < R, Eq. 32) e zona 3 até R_e (Eq. 33, R_e com L_m = 10 H, `lm=`); f = K⁻¹·v'_opt,
   nulo para r ≥ R_e, fora do solo e para h_w = 0. Inclui o ótimo degenerado m → 0 dos taludes íngremes (Φ(m)
   crescente a partir de Φ(0⁺) = A: β > 80,76° para α = 1, 81,10° para α = 2, 88,40° para α = 4, nenhum para
   α ≥ 5), em que a zona 1 é o campo tangencial k_h γw sen β / A e_θ. J*(v'_opt) pela Eq. 40 e J*/(k_h H² γw²).
   Numérica: h2 (Eqs. 37–38) por Dormand–Prince 5(4) adaptativo (rtol 10⁻¹³) do estado escalado
   (h2, w = c h2'/m, ∫w²/c, ∫d h2²), bem condicionado para m → 0; m como no Python (Φ(m) em 81 pontos de
   ln m ∈ [ln 10⁻⁷, ln 30], minimização de Brent e raiz de Brent do ponto fixo m = √(C/D)); h2 e h2' tabelados em
   2049 ângulos (o integrador para em cada um) e interpolados por Hermite cúbico com h2' e h2'' da EDO; a integral
   da zona 3 por Gauss–Kronrod adaptativo; só STL. Construção 2–10 ms; avaliação 40–80 ns por ponto, sem estado
   mutável (seguro entre threads). `Force(x, f)` recebe coordenadas NeoPZ e devolve f com y para cima
   (x_artigo = x, y_artigo = −y, f_y trocado de sinal); `ForcePaper`/`VelocityPaper` usam as do artigo;
   `AnalyticalForceField(ptr)` dá o `slope::ForceField`, usado por `fs water=analytical` (γ_ref = γ').

6. **Análise limite cinemática** (`LimitAnalysis.h`, porte de `scripts/limit_analysis.py`, Seção 4.2, Eqs. 42–58):
   limite superior Γ = min P_mr/(P_γ + P_u) (Eq. 55) sobre mecanismos rotacionais em espiral logarítmica,
   H_crit = Γ H. Internamente nas coordenadas do artigo (y para baixo), como o Python; os campos entram em
   coordenadas NeoPZ (`slope::ForceField`) e são convertidos a cada chamada (ponto (x, y_artigo) → (x, −y_artigo),
   f_y trocado de sinal, u inalterado). Parametrização única (θ1, θ2, s): mecanismo I (B na face,
   s = η ∈ (0, 1]) e II (B no terreno do pé, s = 1 + d/H, d ≤ 10 H); r0, C, L pelas Eqs. 43–46; teste de
   admissibilidade exato do Python (0 < θ1 < θ2 < π, r0 > 0 sem cancelamento, r_h ≤ 1000 H, A no topo à esquerda
   de O, C do lado do ar da face, θ2 < π − β em I, espiral abaixo de O e, em II, de T). P_mr pela Eq. 48; P_γ pela
   forma fechada f1 − f2 − f3 (f3 numa forma válida em β = 90°; f4 = triângulo C–T–B em II) e por quadratura;
   P_u = ∫ f·U para qualquer `ForceField` por quadratura de Gauss composta em coordenadas polares em torno de C
   (raios divididos em O e T, painéis graduados nas pontas, raios cortados nos círculos de descontinuidade do campo
   e quebras angulares nas tangentes e interseções desses círculos — r = R_w, R, R_e e círculos graduados para a
   singularidade r^(m−1) em O no campo analítico), ou, para f = −∇u (campo FE), pela fórmula de contorno exata
   P_u = −∮ u U·n ds (div U = 0 numa rotação rígida): só u na superfície A–O–(T)–B e tan φ ∫u r² dθ na espiral
   (`la::PotentialOf(PoreField)` expõe u = o interpolante P2 de `PoreField::Evaluate`). Otimização como no
   Python: para cada classe (I, II) e semente, PSO de melhor global (coeficientes de Clerc–Kennedy, 40 partículas,
   ≤ 150 iterações) com a quadratura `coarse`, polida por Nelder–Mead limitado (o algoritmo do scipy, 2 reinícios)
   com a `fine`. Os números aleatórios são de `std::mt19937_64` (não os do numpy): as trajetórias do PSO diferem,
   os ótimos polidos coincidem. Diferença deliberada: quando menos de n/2 mecanismos do conjunto inicial (4n
   sorteios admissíveis) têm P_ext > 0, mais conjuntos são sorteados (`pools=25`; `pools=0` = Python, que desiste
   com Γ = ∞ se nenhum tem). As rodadas (classe × semente) são independentes e rodam em threads (`threads=`); o
   resultado não depende do número de threads.

Reuso: `SlopeAnalysis.h` e `SlopeModel.h` incluídos sem alteração; `PoreField`, `SetSeepageForce` e
`FactorOfSafety` de `SlopeDrawdown/main.cpp` adaptados (`SeepageForceField.h`, `FEMStability.h`). Nenhum arquivo
de `SlopeDrawdown`, `SlopeMohrCoulomb` ou da biblioteca foi modificado.

## Compilar e executar

```
cmake -S <neopz> -B <build>      # uma vez, com BUILD_PLASTICITY_MATERIALS=ON e BUILD_PROJECTS=ON
ninja -C <build> SlopeSeepageForces
SlopeSeepageForces <comando> [chave=valor ...]
```

| comando | o que faz |
|---|---|
| `mesh` | malhas hidráulica e de estabilidade, qualidade e consistência (`vtk=1` grava; `sweep=1` testa β = 15…90°, h_w/H = 0…1) |
| `verify` | soluções manufaturadas (linear e quadrática K-harmônica, exatas em P2) e localização de pontos |
| `seepage` | J, J/(k_h H² γw²) e min p para `href=0,1,2` (ordem observada e extrapolação de Richardson); `csv=`, `vtk=` |
| `fig5` | J/(k_h H² γw²) para `alphas=` × `betas=` contra as curvas tracejadas da Fig. 5 (`data/fig5_vector_fill_polygons.csv`) e −J*(v'_opt)/(k_h H² γw²) do campo analítico contra as cheias (`fe=0`: só o analítico); `out=<csv>`: produção na grade de `reproduce_python.py` (β = 15…90° de 7,5 em 7,5°, `href=1`) num CSV retomável |
| `analytical` | campo analítico K⁻¹·v'_opt: R_w, R, R_e, A, m, C, D, F, J* (Eq. 40) por zona e f em pontos de amostra; `pts=` (zona e f em pontos `x,y` do artigo), `csv=` (grade), `bench=<n>` (tempo por ponto), `ref=<arquivo>` (comparação com a referência Python de `scripts/analytical_seepage_reference.py`) |
| `probe` | u e f nos pontos de `pts=<arquivo>` (linhas `x,y` em coordenadas do artigo) |
| `fs` | fator de aumento de gravidade com as forças de percolação (`water=seepage`), com o campo analítico (`water=analytical`, `lm=10`) ou seco (`water=dry`) |
| `check` | autotestes (~20 s, código de saída 1 se algum falhar): convenções, dados de contorno, localização, threads e vetor de carga (ver abaixo); campo analítico (ver *(d)*); análise limite (ver *(e)*); drivers das figuras (ver *(f)*); `fembatch` (ver *(g)*) |
| `la` | análise limite cinemática (`LimitAnalysis.h`): Γ = min P_mr/(P_γ + P_u) sobre os mecanismos I e II, H_crit = Γ H; imprime Γ, o mecanismo (θ1, θ2, η ou d/H; A, B, C), P_mr, P_γ, P_u |
| `labatch` | `la` para cada linha de `cases=<arquivo>` (opções `chave=valor`), uma linha por caso acrescentada a `out=<csv>`; casos já presentes em `out` são pulados (retomável) |
| `fig8` | Fig. 8: H_crit = Γ(H) H × h_w/H (α = 1, H = 1 m), painéis London (β = 30°, 60°) e Israeli (35°, 60°), curvas `vopt` (K⁻¹·v'_opt) e `FE` (−∇u'_FE, caixa 50/10/30 m) e a ponta comum h_w = 0; uma linha por caso em `results/cpp/fig8.csv` (retomável: casos já presentes com as mesmas configurações são pulados), progresso em `fig8.log` |
| `fig9` | Fig. 9: Γ × β (H = 5 m, h_w = H) para α = 1, 5, 10, curvas `vopt` e `FE` (caixa 50/10/30 m); `results/cpp/fig9.csv` e `fig9.log` (retomável) |
| `fembatch` | Γ_FEM (aumento de gravidade, `FEMStability.h`) dos casos da Fig. 9 (`fig=9`: α = 1, 5, 10 × β = 30…90° de 15 em 15° × `FE`, `vopt`) ou da Fig. 8 (`fig=8`: os quatro painéis × h_w/H = 0; 0,2; 0,5; 1), com a análise limite do mesmo caso na linha; `results/cpp/fem_fig9.csv`, `fem_fig8.csv` (retomável; `scripts/run_fem_batch.sh`; ver *(g)*) |

Opções principais (padrão): `H=5 beta=45 hw=1` (h_w/H) `gamma=20 gammaw=9.81 c=10 phi=30 E=20000 nu=0.3`;
`alpha=1 horder=2 hbc=zero_lb` (`impermeable`, `zero_b`, `zero_l`, `zero_lbr`, `toe_r`; ou `hbcleft=`,
`hbcbottom=`, `hbcright=` = `noflow|zero|toe`); malha hidráulica `hmesh=gen|trig hleft=50 hright=10 hdepth=30`
(unidades de H) `hh0=0.025 hhs=0.0625 hgrade=0.15 hhmax=2 href=0`; malha de estabilidade
`smesh=gen|trig sa=2 sleft= sright= sdepth= sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 sref=0`;
`fs`: `form=u|p|p+ nref=3 mark=0.1 srm=0 maxnewton=100 tolfs=0.002 checkforms=1 vtk=<prefixo>`, `water=fe` (= `seepage`)
`|none` (h_w = 0: f = 0 com γ') e `hboxm=` (caixa hidráulica em metros, como em `la`); imprime também λ extrapolado
para h → 0 (ver *(g)*).
`fembatch`: `fig=9|8`; MEF `nref=3 mark=0.05 sa=2` (caixa de estabilidade sa H + H/tan β, limitada pela caixa
hidráulica) `sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 maxnewton=100 tolfs=0.002`; análise limite e campo como em
`fig8`/`fig9` (o mesmo `href=1` para o campo FE do MEF); casos: `fig=9` `H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1
alphas=1,5,10 betas=30,45,60,75,90 curves=vopt,FE`; `fig=8` `soil=swapped gammaw=9.8 panels=… hws=0,0.2,0.5,1
curves=vopt,FE`; `out=<csv>`.
`la`: `water=fe|analytical|none|dry` (campo FE −∇u_h | K⁻¹·v'_opt com `lm=10` | f = 0 com γ' = γ − γw, o caso
h_w = 0 | f = 0 e γw = 0) `hboxm=50,10,30` (caixa hidráulica em metros: à esquerda de O, à direita de T, abaixo de T;
substitui `hleft`…) `href=0 pu=auto|domain|boundary` (auto: contorno para o campo FE) `mech=I,II seeds=0,1,2 np=40
niter=150 pools=25 qsearch=coarse qfinal=fine polish=1 dmax=10 threads=4 verbose=0 out=<csv>` (acrescenta uma
linha) `x=θ1,θ2,s` (sem otimização: P_mr, P_γ e P_u desse mecanismo por todas as regras).
`fig8`/`fig9`: análise limite com as opções de `la` (`seeds=0,1,2 np=40 niter=150 pools=25 qsearch=coarse qfinal=fine
polish=1 dmax=10 mech=I,II lm=10`, mas `threads=2`), caixa hidráulica em metros `hboxm=50,10,30` (`hbc=zero_lb`),
`href=1`, `curves=vopt,FE`, `out=<csv>` (progresso em `<csv>.log`). `fig8`: `soil=swapped|table1` (pares (c, φ) da
Tabela 1 trocados entre os painéis | como impressos) `gammaw=9.8` `H=1` `panels=London30,London60,Israeli35,Israeli60`
`hws=0,0.05,0.1,0.2,…,1` (padrão `results/cpp/fig8.csv`, ou `fig8_table1.csv` com `soil=table1`; c, φ, γ = 18, β e
h_w vêm dos painéis). `fig9`: `H=5 c=10 phi=30 gamma=20 gammaw=9.81 hw=1 alphas=1,5,10 betas=15,20,…,90` (com 37,5,
52,5, 67,5 e 82,5; padrão `results/cpp/fig9.csv`). As linhas guardam as colunas de `data/paper_fig*.csv`, depois
Γ/H_crit, os parâmetros, o mecanismo e uma cadeia `settings`; uma linha cortada por uma interrupção é ignorada e o
caso é refeito.
`hmesh=trig`/`smesh=trig`: `TriGMesh(1 + ref)` de `SlopeMohrCoulomb` (H = 10, β = 45°, caixa 70 × 40 m).

Exemplos:

```
SlopeSeepageForces check
SlopeSeepageForces verify alpha=5
SlopeSeepageForces seepage beta=45 alpha=5 href=0,1,2
SlopeSeepageForces fig5 alphas=1,4 betas=30,60,90 href=1
SlopeSeepageForces fig5 fe=0 betas=15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90   # curvas cheias (analítico)
SlopeSeepageForces analytical beta=30 alpha=1 H=1 bench=1000000
SlopeSeepageForces analytical ref=data/analytical_seepage_reference.csv
SlopeSeepageForces fs beta=60 nref=3                       # dados da Fig. 9, h_w = H, alpha = 1
# regressão contra SlopeDrawdown (talude de SlopeMohrCoulomb, procedimento idêntico):
SlopeSeepageForces fs H=10 beta=45 gammaw=10 water=dry smesh=trig srm=1 nref=3 maxnewton=30
SlopeSeepageForces fs H=10 beta=45 gammaw=10 hbc=impermeable hmesh=trig horder=1 form=p+ smesh=trig srm=1 nref=3 maxnewton=30
# análise limite: Fig. 9 (beta = 60, alpha = 1) com o campo FE na caixa do artigo em metros e com o analítico;
# ponta h_w = 0 da Fig. 8 (painel London, parâmetros trocados); corte vertical seco (N = 3,8313)
SlopeSeepageForces la H=5 beta=60 hw=1 water=fe hboxm=50,10,30
SlopeSeepageForces la H=5 beta=60 hw=1 water=analytical verbose=1
SlopeSeepageForces la H=10 beta=60 hw=0 water=none c=11.7 phi=24.7 gamma=18 gammaw=9.8
SlopeSeepageForces la H=1 beta=90 c=1 phi=0 gamma=1 water=dry
SlopeSeepageForces labatch cases=results/limit_analysis_cpp/cases_vs_python.txt out=minha_saida.csv   # retomável
# produção das Figs. 5, 8 e 9 (retomável, ~25 min com 2 threads) e figuras/tabelas de comparação:
sh scripts/run_cpp_figures.sh <build>/Projects/SlopeSeepageForces/SlopeSeepageForces 2
# ou, passo a passo (a partir deste diretório):
SlopeSeepageForces fig5 out=results/cpp/fig5.csv
SlopeSeepageForces fig9 threads=3                         # results/cpp/fig9.csv
SlopeSeepageForces fig8                                   # soil=swapped gammaw=9.8 H=1 -> results/cpp/fig8.csv
SlopeSeepageForces fig8 soil=table1                       # Tabela 1 como impressa -> results/cpp/fig8_table1.csv
SlopeSeepageForces fig9 alphas=5 betas=60 curves=FE out=teste.csv
# MEF (aumento de gravidade) como verificação da análise limite (seção (g)):
SlopeSeepageForces fs H=5 beta=60 hboxm=50,10,30 href=1 sright=2 nref=3 mark=0.05      # Fig. 9, alpha = 1, campo FE
sh scripts/run_fem_convergence.sh                 # estudo de convergência -> results/fem/convergence/ (retomável)
python3 scripts/fem_convergence_table.py          # tabela do estudo
CPUS=2,3 sh scripts/run_fem_batch.sh "" 9,8       # produção -> results/cpp/fem_fig9.csv, fem_fig8.csv (retomável)
python3 scripts/fem_batch_table.py                # MEF x análise limite x artigo
python3 scripts/python_fig9_box_metres.py                 # referência Python da curva FE da Fig. 9 (caixa em metros)
python3 scripts/plot_results.py                           # results/fig{5,8,9}.png e results/cpp/comparison_*
```

## Resultados (verificação)

Tempos de parede em 4 CPUs compartilhadas (montagem elastoplástica com 4 threads), build Release.

### (a) Percolação

* **Soluções manufaturadas** (`verify`; Dirichlet na superfície, fluxo exato K∇u·n nas laterais e na base, domínio
  50/10/30 H): u linear e u quadrática K-harmônica (c₃(k_v x² − k_h y²) + c₄xy) reproduzidas em arredondamento —
  P2, β = 45°, α = 5 (9165 equações): erro relativo máximo de u 1,8·10⁻¹⁴ e 2,6·10⁻¹⁴, de ∇u 2,4·10⁻¹¹ e
  5,0·10⁻¹³, de J 9,6·10⁻¹⁵ e 1,9·10⁻¹⁴; β = 90°, α = 10, h_w = 0,4 H: u 1,6·10⁻¹⁴; P1 linear, β = 15°, α = 3:
  6,5·10⁻¹⁵ (0,3–0,7 s por caso). Localização de pontos: 10⁶ pontos aleatórios, 0 não encontrados no solo,
  0 encontrados fora, erro de u 1,8·10⁻¹⁴, 0,2 µs por ponto.
* **Malhas** (`mesh sweep=1`): 224 malhas (β = 15…90° de 5 em 5°, h_w/H = 0…1, os dois domínios), área e
  comprimentos de contorno por id exatos a 6·10⁻¹⁵, conformes (cada aresta de triângulo é compartilhada por dois
  triângulos ou por um triângulo e uma linha de contorno), ângulo mínimo 28,0°, 4 s no total.
* **Autotestes** (`check`, 99 verificações, 6 s; sem vazamentos nem erros de memória no valgrind):
  * soluções manufaturadas com α = 5 (e o mesmo problema com K transposto, que tem de falhar: confirma
    K = diag(k_h, k_v) em (x, y));
  * no campo resolvido (β = 15, 30, 45, 60, 90°, h_w/H = 0,25…1, α = 1…10, condições `zero_lb`, `toe_r`,
    `impermeable`, `zero_b`): cada linha de contorno está na sua parte do contorno (ids −1…−6); u nos pontos das
    linhas de topo, face e pé igual à Eq. 21 escrita de forma independente em coordenadas do artigo (erro ≤ 9·10⁻¹³
    γw H) e igual ao dado das laterais com Dirichlet; p ≥ 0; f = −∇u (diferenças centrais de u, 2·10⁻¹⁰);
    adaptador `FromPaperCoordinates`; vértices, pontos médios e centroides de todos os triângulos localizados, com o
    mesmo valor P2 vindo de qualquer triângulo vizinho (continuidade, 6·10⁻¹⁶); nenhum ponto a 10⁻⁶ H fora do solo
    encontrado; avaliador chamado de 4 threads bit a bit igual ao serial;
  * vetor de carga do problema de estabilidade: a resultante (vetor de carga aplicado às translações rígidas) é
    igual a −∮u n ds − γ' A e_y (teorema da divergência) a 8·10⁻¹⁵ com as malhas `trig` (integração exata) e a
    3·10⁻⁶ com as do gerador (erro de quadratura de f descontínuo); F(λ = 2,7) = 2,7 F(1) a 1·10⁻¹⁶ (λ multiplica
    γ' e f); formas u e p iguais a 3·10⁻¹⁷; montagem com 4 threads igual à serial; nenhum ponto de integração da
    malha de estabilidade fora da malha hidráulica (`fs` agora recusa esse caso, em que f seria zero sem aviso).
* **Convergência de J/(k_h H² γw²)**, β = 45°, h_w = H, caixa 50/10/30 H, refinamentos uniformes:

| `hbc` | α | href 0 (9165 eq., 0,5 s) | href 1 (36201 eq., 3,3 s) | href 2 (143889 eq., 36 s) | ordem em h | extrapolado |
|---|---|---|---|---|---|---|
| zero_lb | 1 | 0,7531851705 | 0,7531749864 | 0,7531720833 | 1,81 | 0,75317093 |
| zero_lb | 5 | 0,3103064392 | 0,3102939854 | 0,3102902683 | 1,74 | 0,31028869 |
| impermeable | 1 | 0,6965320306 | 0,6965220322 | 0,6965191815 | 1,81 | 0,69651804 |
| impermeable | 5 | 0,2994644108 | 0,2994521552 | 0,2994484974 | 1,74 | 0,29944694 |

  J converge por cima (espaço conforme com o dado de Dirichlet exato); a ordem < 2 vem da singularidade no pé
  (canto reentrante de 180° + β). Mínimo de p = u − γw y em todos os nós P2 e centroides: −1·10⁻¹¹ kPa
  (arredondamento, na aresta do topo): p ≥ 0.
* **Fig. 5** (curvas tracejadas, valores vetoriais de `data/fig5_vector_fill_polygons.csv`), href = 1:

| α | β | artigo | `zero_lb` | dif. | `impermeable` | dif. |
|---|---|---|---|---|---|---|
| 1 | 30 | 0,68606 | 0,685846 | −0,03 % | 0,627196 | −8,6 % |
| 1 | 60 | 0,80848 | 0,807917 | −0,07 % | 0,752123 | −7,0 % |
| 1 | 90 | 0,91037 | 0,908746 | −0,18 % | 0,853657 | −6,2 % |
| 4 | 30 | 0,32594 | 0,325780 | −0,05 % | 0,312898 | −4,0 % |
| 4 | 60 | 0,36835 | 0,368011 | −0,09 % | 0,355331 | −3,5 % |
| 4 | 90 | 0,39565 | 0,395027 | −0,16 % | 0,382354 | −3,4 % |

  Com `zero_lb` a grade completa da Fig. 5 (α = 1, 2, 4, 10 × β = 15…90° de 15 em 15°, href 0, 12 s) fica entre
  −0,01 % e −0,22 % do artigo (o artigo um pouco acima, como esperado de uma discretização mais grossa);
  em β = 45° (href 0) as outras predefinições ficam, para α = 1 / 10: `zero_b` −0,06 % / −1,0 %, `zero_l`
  −7,3 % / −0,5 %, `toe_r` +19 % / +100 %, `zero_lbr` +174 % / +158 %.
  Logo o domínio do artigo tinha u = 0 à esquerda e na base (as isolinhas da Fig. 4b terminam na lateral direita,
  impermeável). Conferência com a implementação independente `scripts/fe_seepage.py` (P2, outra malha): J
  `zero_lb` 0,7531709 (Python, ref 2) × 0,75317093 (extrapolado aqui); nos pontos de
  `data/fe_seepage_reference_points.csv` (10 casos: as duas condições, α até 10, h_w = H/2) J difere ≤ 6·10⁻⁵,
  u ≤ 5·10⁻⁵ e f ≤ 0,16 % (relativos).

### (b) Regressão contra SlopeDrawdown (talude de SlopeMohrCoulomb: H = 10 m, β = 45°, caixa 70 × 40 m)

γ = 20, γw = 10, c = 10 kPa, φ = 30°, E = 20000 kPa, ν = 0,3; GI e SRM com marcação das duas zonas plásticas
(`srm=1`), `maxnewton=30` (o procedimento de `SlopeDrawdown`). FS no ciclo 3:

| caso | malha de estabilidade (equações no ciclo 3) | campo hidráulico | FS GI | FS SRM | tempo |
|---|---|---|---|---|---|
| seco | `TriGMesh(1)` (1814) | — | **1,90625** | **1,22852** | 2,0 min |
| seco | gerador, 229 triângulos (4832) | — | 1,85938 | 1,21875 | 5,2 min |
| permanente, base e laterais impermeáveis | `TriGMesh(1)` (1710) | P1 em `TriGMesh(1)`, `form=p+` | 0,77612 | 0,88574 | 5,2 min |
| idem | `TriGMesh(1)` (1826) | P2 convergido (gerador, 4284 eq.), `form=u` | 0,74609 | 0,86523 | 4,3 min |
| idem | gerador (4982) | P2 convergido, `form=u` | 0,72070 | 0,85059 | 10,8 min |

* **Seco**: idêntico ao README de `SlopeDrawdown` (1,906 / 1,229 no ciclo 3; GI 3,04077 / 2,51562 / 2,10254 /
  1,90625 nos ciclos 0–3, o ciclo 0 conferido também rodando `SlopeDrawdown nref=0`), com a mesma sequência de
  malhas de `SlopeMohrCoulomb` (870 / 918 / 1140 / 1814 equações). Com o gerador o FS é menor no
  mesmo ciclo porque a malha inicial já é graduada junto à face (mais equações): ambas convergem para os valores
  refinados de `SlopeMohrCoulomb` (ν = 0,3, ciclo 5: 1,789 / 1,207) e de Bishop (1,806 / 1,206).
* **Percolação permanente** (o regime "steady" de `SlopeDrawdown`, p = 0 em toda a superfície = Eq. 21 com
  h_w = H): rodando o próprio `SlopeDrawdown nref=0` (4 min) o ciclo 0 do estado `steady` é GI 1,12439 /
  SRM 1,0625, **idêntico** ao desta implementação com o campo P1 na mesma malha (e o seco 3,04077 / 1,39868).
  No ciclo 3 o README dá 0,781 / 0,886 e aqui 0,776 / 0,886 (SRM igual). O p de `SlopeDrawdown` vem da análise
  u-p no último passo e difere do Laplace P1 por ~5·10⁻⁴ kPa (termo de armazenamento); perto do colapso o driver
  com `maxnewton=30` amplifica perturbações desse tamanho: mudando γw de 10 para 10,00005 (5·10⁻⁶) o GI vai de
  0,83252 para 0,82813 no ciclo 2 e de 0,77612 para 0,77832 no ciclo 3 (SRM inalterado). A diferença de 0,6 % no
  GI é desse ruído do driver, não do campo nem do procedimento. **Confirmado** com um programa de teste (fora do
  projeto) que inclui `SlopeDrawdown/main.cpp` e passa o próprio campo `steady` da análise u-p às duas cadeias, na
  mesma `TriGMesh(1)` com `nref=3`: o `FactorOfSafety` de `SlopeDrawdown` dá GI 1,12439 / 0,936523 / 0,83252 /
  0,78125 e SRM 1,0625 / 0,966309 / 0,912476 / 0,885742 (= README de `SlopeDrawdown`), e `GravityIncreaseFS`
  deste projeto (com `maxnewton=30`) dá exatamente os mesmos números em todos os ciclos (870 / 918 / 1140 / 1710
  equações). O campo u-p difere do Laplace P1 daqui por 1,4·10⁻³ kPa em p (máx. 362 kPa) e 7,6·10⁻⁵ kPa/m em ∇p
  (máx. 15 kPa/m), uma perturbação relativa de 5·10⁻⁶ que muda o GI do ciclo 3 em 0,66 %.
* O campo P1 da malha `TriGMesh(1)` é grosso (J/(k_h H² γw²) = 0,4437 contra 0,4170 convergido, +6 %): com o
  mesmo λ-driver e a mesma malha de estabilidade, o ciclo 0 cai de 1,124 (P1 `TriGMesh(1)`) para 1,104 (P1 href 3)
  e 1,099–1,101 (P2 convergido); no ciclo 3 o campo convergido dá GI −3,9 % e SRM −2,3 %. As forças de percolação
  junto ao pé (singulares) ficam subestimadas no campo grosso.
* **Formas da força de volume** (malhas do gerador, campo P2, `checkforms=1`): os vetores de carga de
  b = λ(γ_sat g − ∇p) e b = λ(γ' g − ∇u) coincidem a |F_u − F_p|/|F_u| = 3,4·10⁻¹⁸ (a parcela de percolação é
  12 % de |F|); nenhum ponto de integração tem p < 0 (min 0,31 kPa), logo p⁺ = p. λ nos ciclos 0 e 1:
  `maxnewton=100` → 0,892578 e 0,795166 nas duas formas, mesma malha refinada (1410 eq.); `maxnewton=30` →
  0,886719 (u) × 0,880859 (p) no ciclo 0: perto do colapso o Newton com retrocesso precisa de 40–100 iterações e
  o limite de 30 declara falha em tentativas que convergiriam, tornando o λ aceito sensível ao arredondamento
  (0,5–1 % abaixo). Por isso o padrão deste projeto é `maxnewton=100` (custo semelhante: menos tentativas
  falhas); `maxnewton=30` reproduz `SlopeMohrCoulomb`/`SlopeDrawdown`.

### (c) Estabilidade com percolação: domínio e tempos

Domínio de estabilidade proposto (`sa=2`): 2 H + H/tan β à esquerda de O, à direita de T e abaixo de T. Para
β = 45° e H = 10 m é exatamente a caixa 70 × 40 m de `SlopeMohrCoulomb`; o termo H/tan β acompanha o
comprimento da face, que governa o tamanho do mecanismo dos taludes abatidos. O domínio hidráulico (50/10/30 H)
contém o de estabilidade para β ≥ 15° (com `sa` até 6; `fs` recusa uma malha de estabilidade com pontos de
integração fora da malha hidráulica, onde f seria zero sem aviso).

Dados da Fig. 9 (H = 5 m, c = 10 kPa, φ = 30°, γ = 20, γw = 9,81, h_w = H, α = 1, `zero_lb`), aumento de gravidade
apenas, `nref=3`, `mark=0.1`, `maxnewton=100` (Γ_FEM = λ; tempos por ciclo, incluindo a percolação, 0,5 s):

| β | caixa (m) | ciclo 0 | ciclo 1 | ciclo 2 | ciclo 3 | total |
|---|---|---|---|---|---|---|
| 30° | x ∈ [−18,7; 27,3], y ≥ −23,7 | 3,16272 (1342 eq., 33 s) | 2,83203 (2006, 44 s) | 2,68262 (3500, 96 s) | 2,59473 (7038, 280 s) | 454 s |
| 60° | x ∈ [−12,9; 15,8], y ≥ −17,9 | 1,16943 (872, 18 s) | 1,09473 (1066, 27 s) | 1,06543 (1354, 24 s) | 1,05371 (1858, 82 s) | 152 s |
| 90° | x ∈ [−10; 10], y ≥ −15 | 0,75195 (682, 18 s) | 0,65527 (914, 20 s) | 0,61230 (1420, 31 s) | 0,58398 (2606, 153 s) | 222 s |

Para referência, a curva FE (−∇u'_FE) digitalizada da Fig. 9 (`data/paper_fig9.csv`) dá Γ = 2,32 / 0,956 / 0,555:
os valores do ciclo 3 ainda estão 12 / 10 / 5 % acima e decrescendo com o refinamento.

* **Sensibilidade à caixa** (β = 30°, ciclo 2): `sa` = 1 / 1,5 / 2 / 3 → λ = 2,71777 / 2,67383 / 2,68262 /
  2,69141; `sa=2` com 6 H abaixo de T: 2,70020. Variação ±0,8 %, não monótona: é ruído da malha inicial (cada caixa
  gera outra triangulação de Delaunay junto ao talude: com a mesma caixa `sa=2` e tamanhos 4 % menores,
  `sh0=0.24 shs=0.24`, λ = 2,65625, −1,0 %), não efeito da caixa. A zona plástica com 1 % do máximo no colapso
  fica, em todas essas caixas, em x ∈ [−0,62 H; x_T + 1,23 H] e até 0,65 H abaixo do pé no ciclo 0 (≤ 0,4 H nos
  seguintes): a pelo menos 1,8 H de qualquer lado mesmo com `sa=1`.
* β = 90° (ciclo 2): `sa` = 2 / 2,5 / 3 / 4 → 0,61230 / 0,60742 / 0,63184 / 0,62158 e, com `sa=2` e tamanhos 4 %
  menores, 0,60547: ±2 %, não monótono, com ruído de malha da mesma ordem (1,1 %); a malha inicial do corte
  vertical é grossa (163–267 triângulos) e no ciclo 0 a dispersão chega a ±6 %. A zona de 1 % vai até 1,0 H atrás
  de O no ciclo 0 e ≤ 0,83 H no ciclo 2 (caixa `sa=2`: 2 H). Em β = 60° vai até 0,89 H atrás de O (folga 1,7 H).
* Conclusão: com `sa=2` o mecanismo fica a ≥ 1 H dos lados (≥ 1,7 H para β ≤ 60°) e caixas maiores não reduzem λ
  além do ruído da malha (1–2 % no ciclo 2). O ruído vem da malha inicial grossa junto ao pé: para a produção
  convém uma malha inicial mais fina (`sh0`, `shs`) ou mais ciclos.
* **Marcação da zona plástica**: com a força de percolação singular no pé, o máximo de √J₂(εᵖ) fica no pé e a
  marcação de 10 % (`SlopeDrawdown`) pode refinar só essa região: em β = 60° a zona marcada no ciclo 3 é
  x ∈ [1,24; 2,89] m, y ∈ [−5,02; −4,07] m e λ estaciona em 1,054. Com `mark=0.02` a marcação cobre toda a
  superfície de ruptura (x ≥ −2,1 m) e λ = 1,16943 / 1,08484 / 1,03687 / 1,00781 (3992 eq. no ciclo 3, 306 s).
  Para a Fig. 9 recomenda-se `mark=0.02` (ou um indicador de dissipação), a decidir na produção (decidido em *(g)*,
  com o campo da caixa em metros: `mark=0.05`, 4 ciclos e extrapolação para h → 0).
* O pivô nulo da skyline (`TPZSkylMatrix::DecomposeLDLt zero pivot`) que aparece às vezes é de uma tentativa muito
  acima do colapso (tangente singular com pontos no ápice); o Newton falha e o passo é reduzido.

### (d) Campo analítico K⁻¹·v'_opt (`AnalyticalSeepage.h`)

* **Contra o Python** (`scripts/analytical_seepage.py`, verificado de forma independente): referência de
  `scripts/analytical_seepage_reference.py` em 13 casos — β = 15, 20, 30, 45, 60, 75, 80, 83, 85, 88 e 90°
  (degenerados: 83°/α = 2, 85°/α = 1, 90°/α = 1; perto do limiar: 80°/α = 1 com m = 0,037 e 88°/α = 4 com
  m = 0,014), α = 1, 2, 4, 5, 10, h_w/H = 0; 0,05; 0,5; 1, H = 1, 5, 10 m, γw = 9,8 e 9,81 —, em pontos aleatórios
  em torno de O (r log-uniforme de 10⁻⁴ H a 1,3 R_e, θ ∈ [0, π]), numa caixa em torno do talude e a 10⁻⁸…10⁻²
  (relativo) dos círculos r = R_w, R, R_e. O f do C++ é avaliado pelo `ForceField` (coordenadas NeoPZ) e convertido
  de volta, de modo que a conversão de coordenadas faz parte da comparação. Com 10⁵ pontos por caso
  (`--npts 100000`, 1,3·10⁶ pontos, 93 MB, não guardado; `analytical ref=<arquivo>`, 1 s; saída em
  `results/analytical_seepage/cpp_port_vs_python_1e5.txt`): nos 838108 pontos com f ≠ 0 a diferença relativa
  máxima |f − f_py|/|f_py| é **8,1·10⁻⁷**, nos 461892 pontos com f_py = 0 (fora do solo, r ≥ R_e, h_w = 0) o C++
  dá f = 0 exato; nenhum ponto excluído (a 10⁻⁹ R_k dos círculos ou a 10⁻¹⁰ H da superfície). m difere em
  ≤ 2,7·10⁻¹⁰, F e J* em ≤ 5·10⁻¹⁴, e a classificação degenerada é a mesma. Por caso, a diferença máxima fica em
  2·10⁻¹³…6·10⁻⁹, exceto nos dois casos de m pequeno (8,1·10⁻⁷ e 4,4·10⁻⁷), sempre no ponto de estagnação da célula
  de recirculação junto ao topo (s* = m^(1/(1−m)), θ ≈ 0, |f| ≈ 0,005 γw), onde f é muito sensível a m. É o erro
  do próprio Python (m pelo LSODA com rtol 10⁻¹²): com DOP853 a rtol 10⁻¹³ o m do Python passa de 0,0372366548918
  a 0,0372366546186 (C++: 0,0372366546186) e, para 88°/α = 4, de 0,0143388911224 a 0,0143388909546 (C++:
  0,0143388909544). O m do C++ muda menos de 10⁻¹⁴ entre rtol 10⁻¹² e 10⁻¹⁴.
* **Fig. 5** (curvas cheias = bordas inferiores das faixas de `data/fig5_vector_fill_polygons.csv`, h_w = H,
  L_m = 10 H; `fig5 fe=0`, 0,4 s para os 64 pontos, `results/analytical_seepage/cpp_port_fig5.txt`): α = 1, 2, 4, 10
  e β = 15…90° de 5 em 5°, diferenças entre −0,020 % e +0,006 % (o Python: 0,02 %). Ótimo degenerado em
  β = 85 e 90° para α = 1 e 2 e em 90° para α = 4.
* **Autotestes** (parte de `check`, 23 verificações, ~1 s): referência Python guardada
  (`data/analytical_seepage_reference.csv`, os 13 casos com 300 pontos cada: 1,6·10⁻⁸; m 2,7·10⁻¹⁰; F e J*
  5·10⁻¹⁴; f = 0 exato onde f_py = 0); os 64 pontos da Fig. 5 (≤ 0,05 %: 0,0198 %); a classificação degenerada igual
  ao critério A > P, P = ∫₀^Θ (∫₀ᵗ d)²/c dt (primeira ordem de Φ(m) = A + m (A − P)), em 86 casos (α = 1, 2, 4, 5,
  10, β = 15…90° e ±0,05° em torno dos limiares 80,76/81,10/88,40°); J* por quadratura do próprio campo pelo
  `ForceField` — ½∫f·K·f no solo até R_e mais ∫u^d (K f)·n na face e no terreno do pé, com a Eq. 21 em coordenadas
  NeoPZ e normais cartesianas, independente da Eq. 40 — igual à Eq. 40 a ≤ 3·10⁻¹¹ (7 casos, α até 10, degenerados
  incluídos); div(K f) = 0 por diferenças centrais (2·10⁻⁸ de |K f|/r); v·e_r contínuo em R_w e R; h2 e h2'
  interpolados contra integrados entre os nós (9·10⁻¹³); condição natural na face c h2' h2 = √C (√C + √D)
  (9·10⁻¹³); f = 0 fora do solo, para r ≥ R_e e com h_w = 0; avaliação de 4 threads igual à serial bit a bit.
  Testes de mutação: trocar o sinal de f_y na conversão para NeoPZ ou pôr α em f_x em vez de f_y (K⁻¹ transposto)
  faz falhar a comparação com o Python, a quadratura de J* (diferenças 1,2–7,6) e o divergente.
* **Custo**: construção 2–10 ms (81–110 integrações de h2), avaliação 36–67 ns por ponto pelo `ForceField`
  (`analytical bench=`; 79 ns numa caixa em torno do talude com β = 30°).
* **Estabilidade**: `fs beta=60 water=analytical nref=0` (dados da Fig. 9, α = 1): λ = 1,70264 no ciclo 0 (872
  equações, 20 s) contra 1,16943 com o campo de elementos finitos na mesma malha — razão 1,46, próxima da razão
  entre as curvas tracejada e cheia da Fig. 9 do artigo em β = 60° (1,4228/0,9557 = 1,49, análise limite).

### (e) Análise limite (`LimitAnalysis.h`)

* **Porte**: em 6 mecanismos fixos (β = 30…90°, φ = 0…32°, I e II) P_mr, P_γ (fechada e por quadratura),
  P_u por quadratura de domínio (potencial polinomial e campos com salto num círculo centrado em O e noutro não
  centrado, nos níveis coarse/fine/xfine, com as divisões nos círculos) e pela fórmula de contorno coincidem com
  `limit_analysis.py` a ≤ 2,4·10⁻¹⁴ (relativo): mesmas regras de quadratura, ponto a ponto. Os 45 números de
  estabilidade secos N = γH_c/c de `results/limit_analysis/dry_stability_numbers.csv` (β = 30…90°, φ = 0…40°,
  incluindo o ótimo do mecanismo II em d = d_max para φ = 0, β ≤ 45° e os três Γ = ∞ de β ≤ φ) são reproduzidos
  a ≤ 8·10⁻⁸ (as 6 casas da tabela), em 0,4 s no total; corte vertical com φ = 0: **N = 3,831337** (Chen 1975:
  3,83).
* **Autotestes** (parte de `check`, 22 verificações, ~1,5 s): 403 mecanismos admissíveis aleatórios (β = 30/60/90°,
  φ = 0/20/32°, I e II, r_h ≤ 20 H): P_γ fechada × quadratura polar `xfine` **2,4·10⁻¹³**, × integração
  independente do polígono A–O–(T)–B + espiral em cordas (fórmulas de *shoelace*, 2001/4001 cordas, Richardson)
  7,3·10⁻¹², Eq. 53 impressa × forma em senos 4,2·10⁻¹², P_u de um potencial polinomial pelo domínio × pelo
  contorno 3,8·10⁻¹², nenhum ponto da espiral acima do terreno; P_mr × quadratura de c r² dθ 2·10⁻¹⁶ e contínuo em
  φ → 0; os valores do Python nos mecanismos fixos (1,7·10⁻¹⁵); 5 números secos (≤ 4·10⁻⁸) e o corte vertical;
  Γ = ∞ para β < φ sem percolação; as pontas h_w = 0 da Fig. 8; invariância de escala; campo uniforme f = (0, g)
  igual a γ' + g sem campo pelos dois métodos de P_u (≤ 1,4·10⁻¹⁵); 4 threads × serial (bit a bit); campo FE da
  Fig. 9 (β = 60°, caixa 50/10/30 m): P_u pelo domínio (`ref`) × contorno 5,9·10⁻⁶, P_u no ótimo do Python × o
  Python 1,5·10⁻⁴ e Γ × o Python 1,0·10⁻⁴.
* **Fig. 8, h_w = 0** (talude submerso, f = 0, γ' = 18 − 9,8, (c, φ) da Tabela 1 trocados entre os painéis):
  H_crit = 156,527 / 17,9495 / 228,829 / 12,9860 m (London 30°, London 60°, Israeli 35°, Israeli 60°) = o Python
  a ≤ 1·10⁻⁷; o artigo (`data/paper_fig8.csv`) dá 156,555 / 17,949 / 229,319 / 12,986 (≤ 0,21 %; 229,319 é uma
  ponta cortada, extrapolada).
* **Invariância de escala**: Γ(H = 5)/Γ(H = 10) − 2 = 0 (seco), −6·10⁻¹² (campo FE com a caixa em unidades de H)
  e 4·10⁻¹² (campo analítico), β = 60°, h_w = 0,6 H, α = 4.
* **C++ × Python** (`scripts/compare_cpp_la.py`, casos de `results/reproduce_python/cache.jsonl`, sementes 0 e 1
  nos dois; `results/limit_analysis_cpp/`: `cases_vs_python.txt`, `cpp_vs_python.csv`,
  `comparison_cpp_vs_python.csv`, `labatch_cases_vs_python.log`; H_crit comparado):
  * campo FE com a caixa do artigo em metros (50/10/30 m, u = 0 à esquerda e na base): 14 casos da Fig. 9 (H = 5 m,
    α = 1, 5, 10, β = 20…90°, 3 com mecanismo II) e 8 da Fig. 8 (parâmetros trocados, h_w/H = 0,1…1, um mecanismo
    de face com η = 0,21; rodados em H = 1 m, a escala do artigo, contra o Python em H = 10 m com a caixa
    50/10/30 H — o mesmo problema por semelhança). Com a malha hidráulica padrão (`href=0`) **21 casos a
    ≤ 4,3·10⁻⁴** (média 2·10⁻⁴); com `href=1` a ≤ 1,6·10⁻⁴ (média 7·10⁻⁵). O C++ fica sempre um pouco acima e a
    diferença cai com o refinamento: é o erro de discretização do campo FE (as duas implementações usam malhas
    diferentes), não da análise limite — no ótimo do Python (β = 50°, α = 10) Γ do C++ é 1,76871 / 1,76822 /
    1,76803 com `href=0/1/2` contra 1,767938 do Python. Exceção: β = 20°, α = 10 (Γ ≈ 38, fora da escala da
    Fig. 9) com 2,4·10⁻³ / 1,1·10⁻³: o ótimo está na fronteira L → 0 da classe (A em O), onde o Nelder–Mead para
    em pontos diferentes conforme a semente (8 sementes no C++: 37,99 a 38,86, dispersão 2 %; o valor do Python
    com 2 sementes, 38,215, também não está convergido) e P_γ + P_u = 76 é a diferença de −976 e 1052, o que
    amplifica 14 vezes as diferenças de P_u.
  * campo analítico K⁻¹·v'_opt (a mesma função nos dois lados): 13 casos (Fig. 9 com α = 1, 5, 10 e β = 15…90°,
    incluindo o caso difícil β = 15°, Γ = 119,08, mecanismo de face η = 0,128; Fig. 8 com h_w/H = 0,05…0,7 e
    mecanismos de face) a **≤ 3,5·10⁻¹¹**. Os casos em que o Python não acha mecanismo com P_ext > 0
    (Fig. 9 com α = 5 e 10, β = 15 e 20°) dão Γ = ∞ também no C++, mesmo com 25 conjuntos iniciais a mais.
* **Métodos de P_u no campo FE** (β = 45°, α = 1, caixa em metros): otimização com P_u pelo contorno Γ = 1,3262534;
  pela quadratura de domínio (descontínua nas arestas dos elementos) 1,3262498 (`fine`), 1,3262432 (`xfine`),
  1,3262519 (`ref`): ≤ 8·10⁻⁶; o padrão para o campo FE é a fórmula de contorno (exata para u_h).
* **Custo** (4 CPUs compartilhadas; 2 classes × 3 sementes, ~27000 avaliações): campo FE (β = 60°) 2,5 s com 1 thread
  e 0,9 s com 3 (mais 0,16 s da percolação); campo analítico 1,6 s / 0,6 s; seco 0,08 s. O Python leva 7–13 s por
  caso FE e 13–380 s por caso analítico com 2 sementes (o C++ 0,4–13 s com 2 threads).

### (f) Reprodução das Figs. 5, 8 e 9 (C++: `fig5 out=`, `fig8`, `fig9`)

Configuração estabelecida (SPEC, achados 1–3): caixa hidráulica do artigo **fixa em metros** — 50 m à esquerda
de O, 10 m à direita de T, 30 m abaixo de T, u = 0 à esquerda e na base, fluxo nulo à direita (`hboxm=50,10,30`,
`hbc=zero_lb`) —, campo FE P2 com `href=1`; campo analítico com L_m = 10 H; análise limite como `la` (mecanismos I
e II, sementes 0, 1, 2, PSO + Nelder–Mead, `pools=25`). Fig. 8 em H = 1 m (H_crit = Γ·1 m; a caixa é então
50/10/30 H, a mesma do Python em H = 10 m por semelhança), α = 1, γ = 18, pares (c, φ) da Tabela 1 trocados entre os
painéis (London c = 11,7 kPa, φ = 24,7°; Israeli c = 6 kPa, φ = 32°; `soil=table1` usa a tabela como impressa);
Fig. 9 em H = 5 m (a caixa é 10/2/6 H). Grades de `scripts/reproduce_python.py` (Fig. 5: β = 15…90° de 7,5 em
7,5°; Fig. 8: h_w/H = 0; 0,05; 0,1; 0,2; …; 1; Fig. 9: β = 15…90° de 5 em 5° e 37,5, 52,5, 67,5, 82,5°). A ponta
h_w = 0 da Fig. 8 é uma única análise sem percolação (f = 0, γ' = γ − γw), gravada nas duas curvas.

Autotestes (parte de `check`, 4 verificações, ~1 s): CSV retomável (chaves, linha cortada ignorada, a seguinte começa
em linha nova), caixa hidráulica fixa em metros (50/10/30 m para H = 1, 5 e 10 m), conjuntos de solo da Fig. 8, e
rodadas mínimas de `fig9` e `fig8` repetidas: a segunda pula todos os casos e deixa os CSVs inalterados, e a ponta
h_w = 0 aparece nas duas curvas com o mesmo valor.

Saídas em `results/cpp/`: `fig5.csv`, `fig8.csv`, `fig8_table1.csv`, `fig9.csv`, `fig8_hw0_gammaw.csv`,
`fig9_gammaw9.8.csv` (colunas de `data/paper_fig*.csv` seguidas de Γ/H_crit, parâmetros, mecanismo — classe,
θ1, θ2, η ou d/H, L/H, A, B, C em coordenadas do artigo, r0, P_mr, P_γ, P_u, dispersão entre sementes, parâmetros
na fronteira da busca —, avaliações, J normalizado do campo, tempos e `settings`), os `.log` de progresso, a
referência Python da curva FE da Fig. 9 com a caixa em metros (`python_fig9_FE_box_m.csv`, abaixo), as comparações
`comparison_fig{5,8,8_table1,9,9_gammaw9.8}.csv` (C++, Python, artigo, diferenças, mecanismos dos dois),
`comparison_summary.txt` e `runtimes.txt`. Figuras no leiaute do artigo (artigo digitalizado em cinza fino, C++ em
traço grosso com os mesmos estilos de linha, Python em círculos vazados): `results/fig5.png`, `results/fig8.png`,
`results/fig9.png` (e `results/fig8_table1.png`), por `scripts/plot_results.py`.

**Tempos** (`scripts/run_cpp_figures.sh`, 2 threads da análise limite, 4 CPUs compartilhadas com outro processo):

| etapa | casos | tempo |
|---|---|---|
| `fig5 out=results/cpp/fig5.csv` (href 1, ~35 000 equações) | 44 × 2 funcionais | 116 s |
| `fig9` | 120 análises limite | 268 s |
| `fig8` | 92 | 405 s |
| `fig8 soil=table1` | 92 | 380 s |
| pontas h_w = 0 (γw = 9,8 e 9,81, os dois conjuntos de solo) | 16 | 1 s |
| `fig9 gammaw=9.8` | 120 | 346 s |
| `scripts/plot_results.py` | — | 3 s |

Por caso (medianas): campo analítico 0,9 s na Fig. 9 e 1,2 s na Fig. 8 (até 16–23 s nos taludes abatidos e nos
mecanismos de face), campo FE
2,5 s na Fig. 9 (1,0 s da percolação) e 5,4 s na Fig. 8 (3,1 s da percolação na caixa 50/10/30 H). O Python leva
7–13 s por caso FE e 13–381 s por caso analítico. A referência Python que faltava (curva FE da Fig. 9 com a caixa em
metros: `data/python_fig9.csv` usa a caixa 50/10/30 H; o diagnóstico D1 de `reproduce_python.py` cobria só os nós de
5° de 20 a 90°) foi completada por `scripts/python_fig9_box_metres.py` (15 casos, 73 s; mesma função e mesmo cache
de `reproduce_python.py`).

**γw** — pontas h_w = 0 da Fig. 8 (talude submerso, f = 0, γ' = 18 − γw), H_crit em m:

| painel | artigo | γw = 9,8 | γw = 9,81 | Tabela 1 como impressa (γw = 9,8) |
|---|---|---|---|---|
| London 30° | 156,555 | 156,527 (−0,018 %) | 156,719 (+0,104 %) | ∞ (β < φ = 32°) |
| London 60° | 17,949 | 17,9495 (+0,003 %) | 17,9714 (+0,125 %) | 12,986 (−27,7 %) |
| Israeli 35° | 229,319 (ponta cortada, extrapolada) | 228,829 (−0,214 %) | 229,108 (−0,092 %) | 68,248 (−70,2 %) |
| Israeli 60° | 12,986 | 12,9860 (+0,000 %) | 13,0019 (+0,122 %) | 17,9495 (+38,2 %) |

Nas três pontas visíveis γw = 9,8 reproduz o artigo a ≤ 1,8·10⁻⁴ (London 60° e Israeli 60° a 3·10⁻⁵) e 9,81 fica
sistematicamente 0,10–0,13 % acima: **padrão de `fig8`: γw = 9,8**. A Fig. 9 não tem ponta h_w = 0; com γw = 9,8 Γ
muda −0,02…+0,31 % para β ≥ 30° (mediana +0,02 %; até +2,6 % nos taludes abatidos fora da escala, β ≤ 25°, onde
P_γ + P_u é diferença de termos grandes) e o desvio médio para o artigo aumenta um pouco (α = 1: FE +0,65 → +0,72 %,
vopt +0,65 → +0,70 %): **padrão de `fig9`: γw = 9,81** (o valor da legenda da Fig. 4b). Nos dois comandos `gammaw=`
muda o valor.

**Fig. 5** (h_w = H, α = 1, 2, 4, 10, 11 ângulos; `results/fig5.png`): C++ × artigo (polilinhas vetoriais)
−J*(v'_opt) entre −0,021 % e +0,009 %, J(u'_FE) entre −0,225 % e −0,014 % (médias por α −0,06 a −0,11 %; o artigo um
pouco acima, como esperado de uma malha mais grossa); C++ × Python: −J* ≤ 3,5·10⁻⁶ (limitado pelas 6 casas do CSV
Python), J_FE ≤ 6·10⁻⁵ (malhas diferentes).

**Fig. 9** (`results/fig9.png`) — Γ do C++ / do artigo (* = ponta cortada do artigo, extrapolada; — = fora da
escala do artigo, Γ > 5):

| β | α = 1: FE | vopt | α = 5: FE | vopt | α = 10: FE | vopt |
|---|---|---|---|---|---|---|
| 25° | 3,390 / 3,364 | 5,897 / 5,583* | 5,915 / 5,987* | 18,406 / — | 7,532 / — | 50,334 / — |
| 30° | 2,333 / 2,322 | 3,579 / 3,490 | 3,786 / 3,756 | 6,848 / 6,871* | 4,694 / 5,176* | 8,807 / — |
| 45° | 1,326 / 1,312 | 1,904 / 1,887 | 1,847 / 1,837 | 2,683 / 2,655 | 2,130 / 2,101 | 2,897 / 2,866 |
| 60° | 0,957 / 0,956 | 1,429 / 1,423 | 1,207 / 1,198 | 1,780 / 1,766 | 1,300 / 1,288 | 1,833 / 1,823 |
| 75° | 0,728 / 0,720 | 1,025 / 1,023 | 0,851 / 0,843 | 1,216 / 1,207 | 0,888 / 0,879 | 1,223 / 1,214 |
| 90° | 0,551 / 0,555 | 0,601 / 0,605 | 0,608 / 0,609 | 0,736 / 0,735 | 0,624 / 0,624 | 0,780 / 0,778 |

| α | curva | pontos visíveis | C++/artigo − 1: média | máx. \|·\| | C++/Python − 1: máx. \|·\| |
|---|---|---|---|---|---|
| 1 | FE | 18 | +0,65 % | 1,13 % | 3,6·10⁻⁵ |
| 1 | vopt | 17 | +0,65 % | 2,56 % | 4,7·10⁻⁶ |
| 5 | FE | 17 | +0,40 % | 1,15 % | 4,4·10⁻⁴ (β = 20°; demais ≤ 1,1·10⁻⁴) |
| 5 | vopt | 16 | +0,73 % | 3,39 % | 3,3·10⁻⁶ |
| 10 | FE | 16 | +0,51 % | 1,66 % | 7,2·10⁻³ (β = 20°; demais ≤ 2,1·10⁻⁴) |
| 10 | vopt | 16 | +0,69 % | 3,10 % | 3,1·10⁻⁶ |

As duas curvas ficam ~0,5–0,7 % acima do artigo em média (máximo 3,4 %, sempre em β = 30–35° na curva vopt, onde
a curva sobe rápido e a digitalização é menos precisa). Γ = ∞ (nenhum mecanismo com P_γ + P_u > 0) em β = 15° no
campo FE e em β = 15–20° no vopt para α = 5 e 10, como no Python. Mecanismos iguais aos do Python em todos os casos:
pé (I com η = 1, B = T) em quase todos; face (η = 0,13 e 0,76) no vopt com α = 1 em β = 15 e 20°; II (B no terreno
do pé) no FE com α = 5 em β = 20–25° e α = 10 em β = 20–30°. A diferença com o Python no FE (2·10⁻⁵…2·10⁻⁴) é a das
duas malhas hidráulicas; em β = 20°, α = 10 (Γ = 37,94 contra 38,22 do Python, fora da escala) o ótimo é o de
fronteira L → 0 descrito em *(e)* (dispersão entre sementes 1 %).

**Fig. 8** (`results/fig8.png`) — C++/artigo − 1 do H_crit (polilinhas interpoladas em escala log; * = ponta cortada):

| painel | curva | h_w/H = 0,05 | 0,1 | 0,2 | 0,3 | 0,5 | 0,7 | 1 | média | C++/Python máx. |
|---|---|---|---|---|---|---|---|---|---|---|
| London 30° | vopt | +2,9 % | +10,7 % | +2,2 % | +5,9 % | +3,6 % | +2,1 % | +0,6 % | +3,2 % | 3·10⁻⁶ |
| London 30° | FE | +13,4 % | +1,7 % | −0,1 % | +0,7 % | −0,5 % | −0,6 % | +0,4 % | +1,2 % | 4·10⁻⁶ |
| London 60° | vopt | +0,1 % | +0,1 % | +0,3 % | +0,0 % | +0,6 % | +0,5 % | +1,0 % | +0,4 % | 5·10⁻⁶ |
| London 60° | FE | +0,0 % | −0,3 % | +0,1 % | −0,4 % | +0,0 % | −0,3 % | +0,8 % | −0,0 % | 2·10⁻⁵ |
| Israeli 35° | vopt | +7,8 %* | +57,4 % | +15,4 % | +6,7 % | +3,9 % | +1,7 % | +1,9 % | +9,6 % | 4·10⁻⁶ |
| Israeli 35° | FE | −5,0 %* | +14,7 % | −20,2 % | −23,6 % | −26,7 % | −22,1 % | −19,0 % | −18,6 % | 5·10⁻⁵ |
| Israeli 60° | vopt | +0,4 % | +0,4 % | +0,3 % | +0,3 % | +0,8 % | +0,8 % | +0,8 % | +0,6 % | 4·10⁻⁶ |
| Israeli 60° | FE | +0,2 % | −0,3 % | +0,3 % | −0,3 % | +0,2 % | +0,3 % | +1,1 % | +0,2 % | 3·10⁻⁵ |

London 60° e Israeli 60° (as duas curvas) ficam a ≤ 1,1 % do artigo em todo h_w; London 30° FE a ≤ 1,7 % fora de
h_w/H = 0,05. As discrepâncias grandes — Israeli 35° FE 19–27 % abaixo para h_w/H ≥ 0,2, vopt de London 30° e
Israeli 35° acima do artigo para h_w/H = 0,05–0,3 (até +57 %) e London 30° FE em 0,05 — são as do Python
(SPEC, achado 4, em aberto): o C++ coincide com o Python a ≤ 5·10⁻⁵ em todos os 92 casos (≤ 1,2·10⁻⁴ com
`soil=table1`), com os mesmos mecanismos (de face, η < 1, no vopt de London 30° para h_w/H = 0,4–0,7, no de
Israeli 35° para 0,2–0,7 e no FE de Israeli 35° para 0,05–0,3; de pé nos demais). Não são, portanto, erros dos
portes, e sim do modelo/dados do artigo (ou da digitalização nas pontas cortadas). Com a Tabela 1 como impressa
(`results/fig8_table1.png`, `results/cpp/comparison_fig8_table1.csv`) os desvios são de −41 a +101 %: confirma a
troca dos pares (c, φ).

### (g) Fator de estabilidade por elementos finitos: convergência e produção (`fembatch`)

Verificação independente das curvas de análise limite das Figs. 8 e 9 pelo aumento de gravidade de
`SlopeAnalysis.h` (Mohr–Coulomb associado, P2, deformação plana; item 4 do *Modelo*), com a mesma força de volume
b = λ(γ' g + f) e o mesmo campo de percolação da análise limite: −∇u'_FE na caixa do artigo em metros (50/10/30 m,
u = 0 à esquerda e na base, `href=1`) ou K⁻¹·v'_opt; Γ_FEM = λ_crit. Domínio de estabilidade: sa H + H/tan β à
esquerda de O e abaixo de T (`sa=2`) e, à direita de T, o mesmo limitado pela caixa hidráulica (em H = 5 m o lado
direito fica a 10 m = 2 H de T: fora da caixa o campo FE não existe e `fs` recusa a malha); malha inicial de
Delaunay com 0,25 H em O, W, T e na face, crescendo 0,25 por H até 1 H; ciclos de refinamento da zona plástica
(`nref`; elementos com √J₂(εᵖ) ≥ `mark` × máx. no colapso divididos em 4, como em `SlopeDrawdown`). Os valores da
seção *(c)* usam o campo da caixa em unidades de H; aqui a caixa é a do artigo, em metros.

**MEF × análise limite.** A análise limite cinemática com mecanismos rotacionais em espiral logarítmica é um limite
superior do fator de colapso exato do mesmo problema (material associado, mesmas forças): Γ_LA ≥ λ*. O MEF
elastoplástico em deslocamentos converge para λ* por cima: numa malha grossa a cinemática discreta não representa a
banda de cisalhamento (não alinhada com os elementos) e o colapso fica alto; cada ciclo divide h por 2 na zona
plástica e λ cai com ordem ≈ 1 em h (diferenças sucessivas na razão 0,4–0,6). O último λ convergido da continuação
fica, em cada malha, até `tolfs` abaixo do colapso daquela malha. Logo λ_FEM(h) ↓ λ* ≤ Γ_LA: os λ dos ciclos podem
ficar acima de Γ_LA, mas o valor extrapolado para h → 0 não deve passar dele, e Γ_LA − λ* é a folga do limite
superior da classe de mecanismos. É o que se observa: no ciclo 4 os λ ainda estão 1,4–2,0 % acima de Γ_LA (β = 60
e 90°, e K⁻¹·v'_opt) ou já 0,5 % abaixo (β = 30°), e os extrapolados ficam 0–2 % abaixo de Γ_LA — 0–1,6 % em
β = 60–90°, onde a espiral é praticamente o mecanismo ótimo, 1,5–2 % em β = 30° —; com K⁻¹·v'_opt as estimativas vão
de −2,0 a +0,5 %, dentro da incerteza da extrapolação. As zonas plásticas (1 % do máximo) no colapso
são as dos mecanismos de pé da análise limite: começam 0,1–0,2 H atrás do ponto A da espiral ótima e saem no pé
(com K⁻¹·v'_opt a zona inteira fica acima do nível do pé).

**Convergência** (`scripts/run_fem_convergence.sh`, logs em `results/fem/convergence/`, tabela e
`summary.csv` por `scripts/fem_convergence_table.py`): dados da Fig. 9 (H = 5 m, c = 10 kPa, φ = 30°, γ = 20,
γw = 9,81, h_w = H, α = 1), λ por ciclo (equações; tempo do ciclo, 2 CPUs compartilhadas), marcação de 2 % e de 5 %;
h → 0: ordem 1 com os ciclos 3–4 (2 λ₄ − λ₃) / ordem observada com os ciclos 2–4 (Richardson):

| β, campo, marcação | ciclo 0 | ciclo 1 | ciclo 2 | ciclo 3 | ciclo 4 | h → 0 | Γ_LA | artigo |
|---|---|---|---|---|---|---|---|---|
| 30°, FE, 2 % | 3,11987 (1240; 28 s) | 2,64746 (2204; 87 s) | 2,44531 (4666; 160 s) | 2,35742 (10592; 544 s) | 2,32227 (26334; 1737 s) | 2,2871 / 2,2988 (p = 1,32) | 2,3328 | 2,3216 |
| 30°, FE, 5 % | 3,11987 (1240) | 2,64746 (2128) | 2,44531 (4050) | 2,35742 (8770; 529 s) | 2,32227 (20988; 2069 s) | 2,2871 / 2,2988 | | |
| 60°, FE, 2 % | 1,19580 (826; 15 s) | 1,08887 (1200; 24 s) | 1,02441 (2268; 82 s) | 0,99219 (4870; 217 s) | 0,97266 (11356; 606 s) | 0,9531 / 0,9426 (p = 0,72) | 0,9571 | 0,9557 |
| 60°, FE, 5 % | 1,19580 (826) | 1,08813 (1162) | 1,02734 (2014) | 0,99219 (3924; 81 s) | 0,97461 (8384; 236 s) | 0,9570 / 0,9570 (p = 1,00) | | |
| 90°, FE, 2 % | 0,74316 (682; 11 s) | 0,64502 (1036; 22 s) | 0,59815 (1750; 49 s) | 0,57422 (3652; 95 s) | 0,56250 (8186; 131 s) | 0,5508 / 0,5513 (p = 1,03) | 0,5513 | 0,5548 |
| 90°, FE, 5 % | 0,74316 (682) | 0,64868 (946) | 0,59912 (1568) | 0,57617 (3168; 97 s) | 0,56250 (6292; 151 s) | 0,5488 / 0,5424 (p = 0,75) | | |
| 60°, K⁻¹·v'_opt, 2 % | 1,81250 (826; 6 s) | 1,61914 (1090; 19 s) | 1,52539 (1854; 67 s) | 1,47852 (4116; 121 s) | 1,44922 (9930; 353 s) | 1,4199 / 1,4004 (p = 0,68) | 1,4289 | 1,4228 |
| 60°, K⁻¹·v'_opt, 5 % | 1,81250 (826) | 1,62354 (1056) | 1,53125 (1754) | 1,47852 (3420; 126 s) | 1,45508 (7266; 326 s) | 1,4316 / 1,4363 (p = 1,17) | | |

Ciclo 5 (β = 60°, FE, 5 %; rodada `b60_m05_n5`): 0,96484 (17228 eq.; 1383 s); extrapolados dos ciclos 3–5 0,9551
(ordem 1) / 0,9526 (p = 0,85), −0,2…−0,5 % de Γ_LA; as estimativas de ordem 1 dos pares de ciclos 2–3, 3–4 e 4–5 são
0,9570 / 0,9570 / 0,9551. Com 2 % o ciclo 5 tem mais de 22571 equações e cada iteração de Newton leva ~8 s (a
fatoração LDLt skyline de `SlopeAnalysis.h` é serial; o ciclo levaria ~1 h): a rodada foi parada e o log
`b60_m02.log` tem os ciclos 0–4.

Estudos em β = 60° (campo FE; λ nos ciclos 0–4 ou os indicados):

* **Marcação**: 1 % → 1,19580 / 1,08813 / 1,02588 / 0,99219 (5304 eq. no ciclo 3); 2 % e 5 % acima; 10 % →
  1,19580 / 1,09802 / 1,03613 / 0,99805 / 0,98047 (5194 eq.). De 1 a 5 % o λ de cada ciclo é o mesmo a 0,3 % (no
  ciclo 3, 0,99219 nos três), com 10–26 % menos equações a 5 % que a 2 %; a 10 % o λ de cada ciclo fica 0,6–1,1 %
  acima (refina menos a banda). Em β = 30° e 90° e com K⁻¹·v'_opt, 5 % dá os λ de 2 % a 0,6 % (β = 30°: idênticos em
  todos os ciclos) com 20–27 % menos equações no ciclo 4.
* **Malha inicial** 0,125 H (`sh0=shs=0.125`, 1502 eq.): 1,10132 / 1,03247 / 0,99512 / 0,97461 (9108 eq.), isto é, o
  ciclo k desta malha ≈ o ciclo k + 1 da de 0,25 H (+0,2…1,1 %); extrapolados 0,9541 (ordem 1) / 0,9496.
* **Refinamento uniforme** sem ciclos (`sref=1`, `sref=2`): 1,08813 (3138 eq.) e 1,02441 (12226 eq.), iguais aos
  ciclos 1 e 2 adaptativos (1200 e 2268 eq.): a marcação refina tudo o que importa, com 2,6–5,4 vezes menos equações;
  em h (uniforme, h ∝ neq^−1/2) a ordem observada é 0,76.
* **Continuação e Newton** (ciclos 0–2): `maxnewton=200` dá os mesmos λ que 100; `maxnewton=30` (o de
  `SlopeMohrCoulomb`/`SlopeDrawdown`) −0,4 / −0,8 / −0,9 % (tentativas que convergiriam são declaradas falhas);
  `tolfs=0.001` +0,2 / 0 / +0,14 %, `tolfs=0.005` −0,4 / −0,3 / −0,3 %: com 0,002 o erro da continuação é ≤ 0,2 %.
* **Domínio** (ciclo 2; base `sa=2` com o lado direito a 2 H de T): β = 30°: `sa=1.5` / 2 / 3 → 2,39258 / 2,44531 /
  2,41895, lado direito a 1,5 H 2,40137; β = 90°: 0,58789 / 0,59815 / 0,59961 e 0,59815. Variação de ±1–2 %, não
  monótona: é o ruído da malha inicial de Delaunay (no ciclo 0, β = 30°: 2,74 / 3,12 / 2,83 e 2,69), não efeito da
  caixa; com `sa=2` a zona de 1 % fica a ≥ 1 H dos lados e do fundo (nas variantes, β = 30°: x ∈ [−0,36 H;
  x_T + 0,55 H] e até 0,28 H abaixo do pé no ciclo 2; β = 90°: até 0,78 H atrás de O). Esse ruído, que cai com os
  ciclos, e a extrapolação limitam a precisão a ~1–2 %.
* Extrapolação: λ × 1/√neq não é linear nestas sequências adaptativas (neq cresce 1,4–2,5 vezes por ciclo enquanto h
  cai pela metade só na zona plástica): a reta pelos dois últimos ciclos dá valores 1–4 % abaixo dos de h ∝ 2⁻ᵏ
  (coluna `lambda_sqrtneq`, só informativa).

**Produção** (padrões de `fembatch`): `nref=3`, `mark=0.05`, `sa=2`, malha inicial 0,25 H, `maxnewton=100`,
`tolfs=0.002`; **Γ_FEM = 2 λ₃ − λ₂** (ordem 1, ciclos 2–3; coluna `Gamma_FEM`, com λ₃ em `Gamma_FEM_last`, todos os
λ dos ciclos, o Richardson dos ciclos 1–3 e o 1/√neq). Contra as extrapolações dos ciclos 3–4 da tabela acima:
β = 30° 2,2695 (−0,8…−1,3 %), 60° 0,9570 (+0,0…+1,5 %), 90° 0,5532 (+0,4…+2,0 %), K⁻¹·v'_opt 60° 1,4258
(−0,7…+1,8 %): precisão de ~1–2 %, enquanto λ₃ sozinho fica 3–4 % acima de Γ_FEM. `nref=4` (Γ_FEM = 2 λ₄ − λ₃)
reduz isso a ~1 % por 2–4 vezes o custo (β = 30°: 48 min por caso). Cada linha guarda também Γ_LA do mesmo caso (`fig9`/`fig8`,
mesmo campo); na Fig. 8 o MEF roda na altura de referência H_ref = H_crit da análise limite (3 algarismos; caixa
50/10/30 H_ref, a do artigo em H = 1 m), para que λ ≈ 1 (a continuação começa com passos de 0,5 e para em 100), e
H_crit = Γ_FEM H_ref — a semelhança λ(H) H = constante é exata para o problema elastoplástico (σ/c, ε e u/H iguais), e
`check` a confirma (vetores de carga em H e 4 H a 6·10⁻¹⁴; λ H de duas rodadas grossas em H = 2 e 8 m a 2 %, a
continuação com `tolfs` 0,02 — 0,07 % com 0,002).

Casos-piloto (`results/cpp/fem_fig9.csv`, `fem_fig8.csv`; 2 CPUs; já contam para a produção):

| caso | Γ_FEM (h → 0) | λ₃ (equações) | Γ_LA | artigo | Γ_FEM/Γ_LA − 1 | tempo |
|---|---|---|---|---|---|---|
| Fig. 9, α = 1, β = 60°, FE | 0,95703 | 0,99219 (3924) | 0,95707 | 0,9557 | −0,004 % | 189 s |
| Fig. 9, α = 1, β = 60°, vopt | 1,42578 | 1,47852 (3420) | 1,42889 | 1,4228 | −0,22 % | 169 s |
| Fig. 8, Israeli 60°, h_w/H = 0,5, FE (H_ref = 3,81 m) | H_crit 3,8193 m | λ 1,0332 (4376) | 3,8111 m | 3,805 m | +0,22 % | 282 s |
| Fig. 8, Israeli 35°, h_w/H = 0,5, FE (H_ref = 8,63 m) | H_crit 8,5879 m | λ 1,01855 (6816) | 8,6303 m | 11,775 m | −0,49 % | 401 s |

No caso Israeli 35° (a curva FE em aberto, SPEC achado 4) o MEF, que não restringe a forma do mecanismo, dá
8,59 m, 0,5 % abaixo da análise limite (8,63 m) e 27 % abaixo do artigo (11,78 m): como a análise limite é um limite
superior do mesmo problema, o valor do artigo não pode ser o colapso destes dados com este campo — a discrepância
não vem da busca de mecanismos.

**Tempo estimado da produção completa** (2 CPUs; por caso, ciclos 0–3 com 5 %: β = 30° 13 min, 60° 3 min, 90°
3 min, mais 1–20 s da análise limite): Fig. 9, 30 casos (≈ 28 min por α e curva) ≈ 3 h; Fig. 8, 28 rodadas (14 nos
painéis de 30/35°, ~13 min, e 14 nos de 60°, ~4–5 min) ≈ 4 h; total ≈ 7 h. Os mecanismos de face da Fig. 8
(K⁻¹·v'_opt de London 30° em h_w/H = 0,5 e de Israeli 35° em 0,2 e 0,5; FE de Israeli 35° em 0,2) são menores que
o talude e podem precisar de mais ciclos.

```
CPUS=2,3 sh scripts/run_fem_batch.sh "" 9,8      # fem_fig9.csv e fem_fig8.csv (retomável: pula os casos feitos)
CPUS=2,3 sh scripts/run_fem_batch.sh "" 9 nref=4 out=results/cpp/fem_fig9_nref4.csv   # ~1 %, 2-4x o custo
python3 scripts/fem_batch_table.py               # MEF x análise limite x artigo -> results/cpp/comparison_fem_fig*.csv
```

Autotestes (parte de `check`, 5 verificações, ~10 s): extrapolações exatas em sequências modelo (ordem 1 e 1,5;
h ∝ neq^−1/2); caixa de estabilidade dentro da hidráulica e igual à regra (sa H + H/tan β limitada) em todos os
casos de produção; semelhança (acima); rodadas mínimas de `fembatch fig=9` e `fig=8` repetidas (a segunda pula todos
os casos; a ponta h_w = 0 nas duas curvas).

## Extensões (a portar)

* Campo analítico K⁻¹·v'_opt: feito (`AnalyticalSeepage.h`, `fs water=analytical`). Para a análise limite, as
  descontinuidades do campo são os círculos r = R_w, R, R_e em torno de O (`Rw()`, `R()`, `Re()`) e, para
  0 < m < 1, |f| ~ r^(m−1) em O (`M()`).
* Análise limite: feita (`LimitAnalysis.h`, comandos `la` e `labatch`; item 6 do *Modelo* e *(e)*).
* Varreduras das Figs. 8 e 9: feitas (`fig8`, `fig9`; *(f)*). Em aberto, como no Python: a curva FE de Israeli 35° e
  as curvas vopt de London 30°/Israeli 35° para h_w/H pequeno (Fig. 8).
* MEF (aumento de gravidade) como verificação independente das curvas das Figs. 8 e 9: `fembatch` (*(g)*), com o
  estudo de convergência e quatro casos-piloto feitos; a produção completa (~7 h com 2 CPUs) é lançada à parte
  (`scripts/run_fem_batch.sh`).

## Limitações

* `hbc=zero_lbr`: no canto superior direito u = 0 (lateral) e u = −γw h_w (pé) são impostos por penalidade nos
  dois elementos de contorno e o nó do canto fica com a média (em `fe_seepage.py` prevalece a superfície).
* `hmesh=trig` não tem nó no nível d'água para 0 < h_w < H (só para a regressão com h_w = H).
* O fator de aumento de gravidade é o último λ convergido: incerteza da continuação ~0,2 % (`tolfs`) mais o ruído
  do Newton perto do colapso (acima).
* O MEF usa a fatoração LDLt skyline serial de `SlopeAnalysis.h` (sem MKL/METIS nesta compilação): acima de ~2·10⁴
  equações cada iteração de Newton leva vários segundos e um ciclo ~1 h, o que limita `fembatch` a 4 ciclos e a
  precisão a ~1–2 % (*(g)*).
