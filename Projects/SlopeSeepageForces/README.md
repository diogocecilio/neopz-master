# SlopeSeepageForces

Estabilidade de taludes sob as forças de percolação de um rebaixamento rápido — reprodução da Seção 4.3 de
Ceron, Cecílio, Linn & Maghous, *Stability analysis of slope subjected to seepage forces considering spatial
variability of soil properties*, IJNAMG 49(11), 2025 (doi:10.1002/nag.3993): Fig. 5 (funcionais hidráulicos),
Fig. 8 (altura crítica × h_w/H) e Fig. 9 (fator de estabilidade × β para α = 1, 5, 10).

Este diretório contém a parte C++ (NeoPZ), construída sobre `Projects/SlopeDrawdown` e
`Projects/SlopeMohrCoulomb`; as implementações de referência em Python estão em `scripts/` (dados digitalizados
em `data/`, resultados em `results/`). A análise limite e o campo analítico K⁻¹·v'_opt
(`scripts/limit_analysis.py`, `scripts/analytical_seepage.py`) ainda serão portados: entram pela interface
`slope::ForceField` (ver *Extensões*).

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
   (`nref`, `mark`; `srm=1` calcula também o SRM e marca as duas zonas, como lá).

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
| `fig5` | J/(k_h H² γw²) para `alphas=` × `betas=` contra as curvas tracejadas da Fig. 5 (`data/fig5_vector_fill_polygons.csv`) |
| `probe` | u e f nos pontos de `pts=<arquivo>` (linhas `x,y` em coordenadas do artigo) |
| `fs` | fator de aumento de gravidade com as forças de percolação (`water=seepage`) ou seco (`water=dry`) |
| `check` | autotestes (~6 s, código de saída 1 se algum falhar): convenções, dados de contorno, localização, threads e vetor de carga (ver abaixo) |

Opções principais (padrão): `H=5 beta=45 hw=1` (h_w/H) `gamma=20 gammaw=9.81 c=10 phi=30 E=20000 nu=0.3`;
`alpha=1 horder=2 hbc=zero_lb` (`impermeable`, `zero_b`, `zero_l`, `zero_lbr`, `toe_r`; ou `hbcleft=`,
`hbcbottom=`, `hbcright=` = `noflow|zero|toe`); malha hidráulica `hmesh=gen|trig hleft=50 hright=10 hdepth=30`
(unidades de H) `hh0=0.025 hhs=0.0625 hgrade=0.15 hhmax=2 href=0`; malha de estabilidade
`smesh=gen|trig sa=2 sleft= sright= sdepth= sh0=0.25 shs=0.25 sgrade=0.25 shmax=1 sref=0`;
`fs`: `form=u|p|p+ nref=3 mark=0.1 srm=0 maxnewton=100 tolfs=0.002 checkforms=1 vtk=<prefixo>`.
`hmesh=trig`/`smesh=trig`: `TriGMesh(1 + ref)` de `SlopeMohrCoulomb` (H = 10, β = 45°, caixa 70 × 40 m).

Exemplos:

```
SlopeSeepageForces check
SlopeSeepageForces verify alpha=5
SlopeSeepageForces seepage beta=45 alpha=5 href=0,1,2
SlopeSeepageForces fig5 alphas=1,4 betas=30,60,90 href=1
SlopeSeepageForces fs beta=60 nref=3                       # dados da Fig. 9, h_w = H, alpha = 1
# regressão contra SlopeDrawdown (talude de SlopeMohrCoulomb, procedimento idêntico):
SlopeSeepageForces fs H=10 beta=45 gammaw=10 water=dry smesh=trig srm=1 nref=3 maxnewton=30
SlopeSeepageForces fs H=10 beta=45 gammaw=10 hbc=impermeable hmesh=trig horder=1 form=p+ smesh=trig srm=1 nref=3 maxnewton=30
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
  Para a Fig. 9 recomenda-se `mark=0.02` (ou um indicador de dissipação), a decidir na produção.
* O pivô nulo da skyline (`TPZSkylMatrix::DecomposeLDLt zero pivot`) que aparece às vezes é de uma tentativa muito
  acima do colapso (tangente singular com pontos no ápice); o Newton falha e o passo é reduzido.

## Extensões (a portar)

* Campo analítico K⁻¹·v'_opt (`scripts/analytical_seepage.py`): implementar um `slope::ForceField`
  (diretamente em coordenadas NeoPZ ou com `FromPaperCoordinates`) e passá-lo a `GravityIncreaseFS` com
  γ_ref = γ'; em `main.cpp`, uma nova opção de `fs` ao lado de `water=seepage|dry`.
* Análise limite (`scripts/limit_analysis.py`): usa `PoreField::Locate`/`EvaluateTri` (ou o `ForceField`) para
  integrar −∇u no mecanismo; `SlopeGeometry` fornece a geometria e `NoSeepage()` o caso h_w = 0.
* Varreduras das Figs. 8 e 9: `Problem::SetSlope(beta, hw/H)` em `main.cpp` reconstrói os dois domínios.

## Limitações

* `hbc=zero_lbr`: no canto superior direito u = 0 (lateral) e u = −γw h_w (pé) são impostos por penalidade nos
  dois elementos de contorno e o nó do canto fica com a média (em `fe_seepage.py` prevalece a superfície).
* `hmesh=trig` não tem nó no nível d'água para 0 < h_w < H (só para a regressão com h_w = H).
* O fator de aumento de gravidade é o último λ convergido: incerteza da continuação ~0,2 % (`tolfs`) mais o ruído
  do Newton perto do colapso (acima).
