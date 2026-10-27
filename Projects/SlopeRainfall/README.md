# SlopeRainfall

Talude homogêneo **não saturado sob chuva** com o modelo acoplado de três fases (esqueleto, água e ar) da
formulação *"Formulação poroelástica não saturada: esqueleto deformável, água e ar"* implementado no NeoPZ, e
reprodução dos exemplos de

> H. Rahardjo, T. H. Ong, R. B. Rezaur, E. C. Leong (2007). *Factors controlling instability of homogeneous soil
> slopes under rainfall*. J. Geotech. Geoenviron. Eng. 133(12), 1532–1543. doi:10.1061/(ASCE)1090-0241(2007)133:12(1532)

É a versão C++/NeoPZ do programa Python de referência `talude3f` (`nsat3.py`, `caso.py`, `rodar_artigo.py`):
mesma geometria, malha, leis constitutivas, condições de contorno, passos de tempo e cálculo do FS.

## Material: `TPZMatUnsaturatedUPwPa` (Material/Plasticity)

Incógnitas (valores totais na solução): deslocamento **u** (H1, P2), pressão da água **p_w** e pressão do ar
**p_a** (H1, P1, manométricas, positivas em compressão; Taylor–Hood). Tração positiva; kN, m, s (kPa).

| Formulação (notas) | Implementação |
|---|---|
| §3 sucção, saturação, p_B | `TPZUnsaturatedHydraulics`: Fredlund–Xing (C(ψ) = 1, como no artigo) ou van Genuchten; Γ = ∂S_w/∂p_w |
| §4 equilíbrio, tensão de Bishop | σ = σ' − α p_B I, p_B = χ p_w + (1 − χ) p_a, χ = S_w (notas) ou constante (χ = tan φᵇ / tan φ' = 1 no artigo); σ' de `TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>` (Mohr–Coulomb do artigo do RHW) |
| §5 Darcy das duas fases | λ_w = **K** k_rw/γ_w, λ_a = (**K**/γ_w)(μ_w/μ_a) k_ra, **K** = R(θ) diag(k₁, k₂) R(θ)ᵀ (anisotropia com eixos principais girados de θ; `SetAnisotropicConductivity(k1, k2, theta)`; isotrópico k₁ = k₂ = k_s); k_rw = S_eᵖ (Leong–Rahardjo) ou Mualem; k_ra de Brooks–Corey ou Mualem; mínimos k_rw,min e k_ra,min (§13.9) |
| §6 fluidos | água com K_w; ar ideal isotérmico, c_a = 1/p_a^abs, termo g_a = c_a q_a·∇p_a mantido |
| §7–8 armazenamento e Biot | S_π[α Δε_v + η Δp_B] + n S_w/K_w Δp_w + ⟨n ΔS_π⟩ + ⟨n S_a/p_a^abs Δp_a⟩, η = (α − n)/K_s |
| §9 porosidade | n constante nos coeficientes (fechamento de 1ª ordem com α = 1 e grãos incompressíveis) |
| §13.6–13.7 Euler implícito e Newton | resíduos das equações de massa multiplicados por Δt; **jacobiano consistente completo** (inclusive ∂λ/∂p, ∂S/∂p, ∂c_a/∂p_a, ∂χ/∂p): convergência quadrática verificada |
| §13.9 estabilidade | Taylor–Hood (LBB) e armazenamento concentrado (⟨·⟩: ∫N_i dΩ vezes o termo com as pressões do vértice, forma secante conservativa n(Sⁿ⁺¹ − Sⁿ)) |

Opções: `SetFrozenPressures` (só o equilíbrio, para o FS com as pressões de um instante), `SetGravityFactor`
(fator do peso e de p_B), `SetGravityScalesSuction`, `SetUnitWeight`/`SetUnitWeightMixture`
(γ = (1 − n)γ_s + n S_w γ_w), `SetLumpedStorage`. As condições de Dirichlet são impostas pela análise (filtro de
equações); `ContributeBC` monta trações e fluxos prescritos de água/ar.

## Problema (Fig. 2 do artigo)

* Geometria: base em y = 0, pé em y = 2 H_s, crista em y = 3 H_s, platôs de 3 H_s; malhas de triângulos do gmsh
  (`scripts/malha_talude.py`, as mesmas da malha "media" do Python: 0,75 m até 2 m de profundidade, crescendo até
  0,3 H_s), guardadas em `meshes/`.
* Solos (Tabela 1): Fredlund–Xing n = m = 1, θ_s = 0,45, k_w = k_s Θ⁴; f10: a = 10 kPa, k_s = 10⁻⁴ m/s;
  f50: 50 kPa, 10⁻⁵ m/s; f100: 100 kPa, 10⁻⁶ m/s. c' = 10 kPa, φ' = φᵇ = 26° (χ = 1), γ = 20 kN/m³ constante.
* Hipóteses de `talude3f` que não estão no artigo: E = 20 MPa, ν = 0,3, ψ = φ'; ar com Brooks–Corey (λ = 3,
  k_ra,min = 10⁻⁴), μ_w/μ_a = 55,6, K_w = 2,2 GPa.
* Contornos: base fixa, laterais em rolete; carga total imposta nas laterais abaixo do NA; base e laterais acima
  do NA impermeáveis; p_a = 0 na superfície (crista, face e pé). NA inicial inclinado de 7° a H_w abaixo do pé,
  pressão hidrostática com sucção limitada a 75 kPa, p_a = 0.
* Chuva (I_r durante 24 h) na projeção horizontal da superfície, condição **unilateral** (conjunto ativo como
  em `unilateral_step`): p_w ≤ 0, infiltração ≤ chuva e complementaridade; vértices com p_w > 0 passam a
  p_w = 0 (escoamento superficial) e voltam a fluxo se a reação pede mais água que a chuva. Sem chuva a mesma
  condição é uma face de percolação.
* Passos: 15 min na chuva, 1 h nas 24 h seguintes, 4 h depois; bissecção do passo (2 níveis, 5 após a ruptura).
* Estado inicial: peso próprio em 4 incrementos com as pressões congeladas; u = 0 mantendo σ'.

### Fator de segurança

Em cada instante de `tfs`, análise independente com p_w e p_a do instante congeladas (como o artigo faz com
SLOPE/W), com a memória do transiente guardada e restaurada:

* **Redução de resistência (SRM)**: estado geostático com a resistência aumentada (F0 = 0,6, para achar FS < 1);
  depois c'/F, tan φ'/F e tan ψ/F aumentam até a perda de convergência (busca de colapso com bissecção até 0,01,
  `collapse_search` do Python). É o FS comparado com o Bishop simplificado do artigo.
* **Aumento da gravidade (GI)**: a partir do estado sem peso, λ multiplica o peso do solo e as pressões positivas
  (peso da água); a sucção é mantida (não é gravitacional; multiplicá-la por λ faria a resistência crescer com o
  peso e o talude nunca romperia).

Se o transiente acoplado perde o equilíbrio com a resistência real (FS < 1), o instante é gravado como ruptura
(`t_fail`) e o transiente continua com a resistência aumentada (F0) só para acompanhar as pressões.

## Compilar e executar

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL=ON -DCMAKE_BUILD_TYPE=Release <neopz>
ninja SlopeRainfall
./SlopeRainfall solo=f50 alpha=45 hs=10 hw=5 I=36 out=f50_I36 threads=4
./SlopeRainfall I=36 passivo=1        # p_a = 0 em todo o domínio (equivale ao SEEP/W do artigo)
./SlopeRainfall I=80 gi=0 vtk=1       # só SRM; VTK (p_w, p_a, sucção, S_w, p_B, u) nos instantes do FS
```
Sem MKL a análise usa `TPZSkylineNSymStructMatrix` (mais lento). Argumentos: `solo`, `alpha`, `hs`, `hw`, `I`
(mm/h), `kx`, `ky` (condutividades saturadas principais, m/s, padrão k_s do solo), `theta` (direção de `kx`, graus anti-horários a partir do eixo x), `train` (h, 24), `ttotal` (h, 48), `tfs` (lista de instantes do FS, h), `passivo`, `gi`, `fs`, `vtk`,
`mesh`, `out`, `verbose` (0–3), `threads`.

Saídas em `out/`: `fs.csv` (t, FS SRM, FS GI), `hist.csv` (chuva, infiltração, vértices alagados, p_w e p_a
máximos, iterações por passo), `snaps_pw.csv` e `snaps_pa.csv` (p_w e p_a nos vértices a cada hora),
`vertices.csv`, `triangles.csv`, `outline.csv`, `summary.txt`.

Séries do artigo (em paralelo) e figuras de comparação:
```
cd scripts
python rodar_series.py --serie B --jobs 4 --exe <build>/Projects/SlopeRainfall/SlopeRainfall
python rodar_series.py --serie todas --jobs 4 --sem-gi
python figuras.py --serie todas            # figuras/<série>_*.png e <série>_tabela.csv
```

## Pós-processamento: FS em log(t) e animação

`scripts/animacao.py <pasta do caso>` (matplotlib) ou `mma/animacao.wl` (Mathematica, adaptado do
pós-processamento do `SlopeDrawdown`: `wolframscript -file animacao.wl <pasta do caso>`) geram

* `fs_tempo.png`: FS por redução de resistência e por aumento da gravidade em escala log(t), t = 0 desenhado em
  0,01 h, linhas em t = 1 h e FS = 1, linha vermelha "ruptura" no instante em que o transiente perdeu o
  equilíbrio;
* `chuva.mp4` e `chuva.gif` (4 quadros/s): p_w preenchido a cada 10 kPa, isolinhas rotuladas a cada 50 kPa,
  p_w = 0 em magenta, setas de chuva enquanto chove; ao lado, o gráfico do FS com o marcador no instante, o FS do
  marcador interpolado linearmente em log(t) (`InterpFS`, constante se houver um só ponto).

Os casos `resultados/B/f50_a45_H10_Hw5_I36` (modelo completo, com FS por aumento da gravidade) e
`resultados/B_passivo_desacoplado/f50_a45_H10_Hw5_I90` trazem os arquivos completos e as animações.

## Resultados

126 casos (malha "media", 4 processos; 15–30 min por caso de 48 h com H_s = 10 m, ≈ 4 h com H_s = 40 m).
Figuras e tabelas em `figuras/` (`python figuras.py --serie todas`); FS, histórico e resumo de cada caso em
`resultados/<série>/<caso>/`. FS = redução de resistência (o do artigo é Bishop simplificado).

Três variantes do modelo na Série B:

* **B** — modelo completo: acoplamento hidromecânico (Biot) e ar ativo;
* **B_passivo** — acoplado, p_a = 0 em todo o domínio;
* **B_passivo_desacoplado** — p_a = 0 e percolação sem o termo S_π α Δε_v (`acoplado=0`): as hipóteses do artigo
  (SEEP/W + análise de estabilidade com as pressões do instante).

### FS inicial

Coincide com o artigo em todas as geometrias: 2,119 (H_s = 10 m, α = 45°; artigo 2,11), 2,547 / 2,350 / 1,775
(α = 26,6 / 33,7 / 63,4°; Tabela 7: 2,54 / 2,38 / 1,69), 3,211 / 1,488 / 1,125 (H_s = 5 / 20 / 40 m; Fig. 10:
3,22 / 1,50 / 1,12), 1,800 … 2,397 com H_w (Fig. 8: 1,80 … 2,35).

### Série B: FS mínimo × I_r (Fig. 7, Eq. 6)

| solo | I_r (mm/h) | Eq. 6 | B (completo) | B_passivo | B_passivo_desacoplado |
|---|---|---|---|---|---|
| f10 | 9 / 36 / 90 | 1,89 / 1,48 / 1,05 | 1,93 / 1,57 / 1,15 | 1,93 / 1,58 / 1,16 | 1,93 / 1,58 / 1,15 |
| f10 | 180 / 360 / 900 | 0,85 / 0,81 / 0,81 | 1,00 / 1,00 / 0,99 | 0,99 / 1,00 / 0,99 | **0,86** / 1,05* / 0,98* |
| f50 | 3,6 / 9 / 18 | 1,98 / 1,79 / 1,49 | 1,96 / 1,74 / 1,32 | 2,01 / 1,86 / 1,62 | 2,01 / 1,85 / 1,61 |
| f50 | 36 / 54 / 90 / 180 | 1,05 / 0,89 / 0,82 / 0,82 | 1,01 / 1,00 / 1,00 / 1,00 | 1,21 / 1,03 / 1,00 / 0,90 | 1,18 / **0,97 / 0,85 / 0,84** |
| f100 | 0,9 / 3,6 / 9 | 2,04 / 1,88 / 1,55 | 2,01 / 1,77 / 1,36 | 2,06 / 1,95 / 1,70 | 2,06 / 1,94 / 1,69 |
| f100 | 18 / 36 | 1,45 / 1,45 | 1,28 / 1,28 | 1,52 / 1,49 | 1,50 / 1,48 |

\* a percolação parou de convergir às 8,5 h com toda a superfície alagada (chuva de 1 e 2,5 k_s no solo mais
permeável): o mínimo está subestimado.

* **Com as hipóteses do artigo** (desacoplado, ar passivo) o FS segue a sigmoide da Eq. 6 nos três solos, inclusive o
  patamar inferior (0,84–0,86 contra 0,81–0,82) e a recuperação após a chuva da Fig. 5b.
* **No modelo completo o FS mínimo não desce abaixo de ≈ 1,00.** O transiente acoplado nunca perde o equilíbrio: com
  Mohr–Coulomb associativo (ψ = φ') a plastificação dilata o solo e, na percolação acoplada, isso gera sucção
  (p_w cai) e aumenta a resistência — o talude se estabiliza no limite. Com ψ < φ' o efeito diminuiria.
* **O ar ativo reduz o FS** nos solos menos permeáveis (f100 com 18 mm/h: 1,28 contra 1,50 com p_a = 0): o ar
  aprisionado sob a frente de umedecimento chega a 50 kPa e eleva p_w (positivo sob a crista, com sucção ainda
  positiva), e a tensão de Bishop com χ = 1 depende de p_w.

### Séries A, C e D (modelo completo)

* **A (Fig. 9)**: FS mínimo sobre os pontos da Fig. 9 nas chuvas fracas e moderadas (f10 9 mm/h: 2,37 / 2,17 / 1,93
  / 1,61 contra 2,35 / 2,16 / 1,88 / 1,45; f50 36 mm/h α = 63,4°: 0,74 contra 0,76; f100 3,6 mm/h); nas chuvas
  fortes, o patamar ≈ 1,0 acima. A Tabela 7 do artigo não reproduz a Fig. 9 em algumas linhas (f10 80 e 360 mm/h,
  f100 80 mm/h); a figura mostra as duas.
* **C (Fig. 8)**: FS inicial e FS mínimo com 9 mm/h sobre as curvas (H_w = 2,5…15 m: 1,38 / 1,74 / 2,01 / 2,17 /
  2,31 contra 1,43 / 1,79 / 2,00 / 2,13 / 2,20); com 36 e 80 mm/h, 1,00–1,14 contra 0,81–1,29.
* **D (Fig. 10)**: H_s = 5 / 10 / 20 / 40 m com 9 mm/h: 2,38 / 1,74 / 1,33 / 1,07 contra 2,48 / 1,79 / 1,37 / 1,08;
  com 36 / 80 mm/h, H_s = 20 m dá 0,86 / 0,99 (artigo 0,97 / 0,88) e H_s = 40 m 1,00 / 0,995 (0,90 / 0,83).

### Anisotropia da permeabilidade

`scripts/anisotropia.py` (figuras `anisotropia_fs.png` e `rotacao_fs.png`): solo f50, α = 45°, 36 mm/h (= k_s),
modelo completo; FS mínimo isotrópico 1,005.

| k_x (m/s) | k_y (m/s) | FS mín SRM | FS mín GI | infiltração em 24 h (m³/m) |
|---|---|---|---|---|
| 1e-5 | 1e-5 | 1,005 | 1,01 | 38,4 |
| 3e-5 | 1e-5 | 1,035 | 1,10 | 40,2 |
| 1e-4 | 1e-5 | 1,078 | 1,25 | 42,7 |
| 1e-5 | 1e-6 | 1,259 | 1,75 | 16,0 |
| 1e-6 | 1e-5 | 1,002 | 1,01 | 35,5 |
| 1e-5 | 1e-4 | 1,195 | 1,56 | 53,9 |

Tensor girado (k_x = 1e-4 na direção θ, k_y = 1e-5 normal a ela; a face desce a −45°, θ = 135° é a direção da face):

| θ (°) | 0 | 15 | 30 | 45 | 60 | 75 | 90 | 105 | 120 | 135 | 150 | 165 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FS mín SRM | 1,08 | 1,19 | 1,33 | 1,33 | 1,28 | 1,24 | 1,20 | 1,16 | 1,15 | 1,18 | 1,20 | 1,11 |
| FS mín GI | 1,24 | 1,64 | 2,21 | 2,16 | 1,95 | 1,73 | 1,55 | 1,42 | 1,37 | 1,49 | 1,60 | 1,36 |
| infiltração (m³/m) | 43,0 | 48,1 | 52,9 | 53,8 | 53,7 | 53,7 | 53,9 | 54,8 | 56,3 | 57,9 | 57,3 | 50,0 |

* k_y < k_s: a chuva infiltra menos (escoa na superfície) e o FS sobe; k_x > k_y espalha a água lateralmente e
  a drena pela face.
* O pior caso girado é θ = 0 (estratificação horizontal, menor infiltração vertical); o melhor, θ ≈ 30–45°
  (direção mais permeável subindo para dentro do maciço, contrária à face): a água se afasta da face.

## Limitações

* Patamar FS ≈ 1 do modelo acoplado com ψ = φ' (acima); a comparação com o artigo deve usar `acoplado=0 passivo=1`.
* Em poucos casos um passo acoplado não converge nem com 5 bissecções (FS ≈ 1,8–1,9): o instante é gravado em
  `t_fail` e as figuras o marcam como "transiente sem convergência" (laranja), não como ruptura.
* Porosidade constante nos coeficientes; sem histerese de retenção; FS por elementos finitos (SRM) e não por
  equilíbrio-limite.
* Custo: a montagem do NeoPZ (elemento multifísico de 3 espaços) domina; Pardiso (MKL) recomendado.

## Alterações na biblioteca

* `TPZMatUnsaturatedUPwPa`, `TPZUnsaturatedHydraulics` (novos, Material/Plasticity).
* `TPZSpStructMatrix`: blocos parcialmente removidos pelo filtro de equações contam só as equações ativas (antes a
  montagem abortava com um componente de um nó prescrito).
* `TPZPardisoSolver`: matrizes não simétricas sempre fatoradas (a iteração CGS com o LU anterior falhava nas
  tangentes de Newton e, com K = 1, também o fallback de fatoração).
