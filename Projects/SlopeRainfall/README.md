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
* Passos: 15 min na chuva, 1 h nas 24 h seguintes, 4 h depois; bissecção do passo (2 níveis; um passo que falha é
  refeito com 5 a partir do início do passo, e depois da ruptura usam-se sempre 5).
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
(`t_fail`) e o transiente continua com a resistência aumentada (F0) só para acompanhar as pressões. Um passo que falha
com 2 níveis de bissecção é tentado de novo com 5, sempre a partir do estado do início do passo (solução, memória e
conjuntos alagado/seco); no modelo desacoplado uma falha com 5 níveis encerra a percolação.

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
`vertices.csv`, `triangles.csv`, `outline.csv`, `summary.txt` (parâmetros, FS inicial e mínimo, `t_fail`, iterações,
`retries`: tentativas de passo descartadas, `retries_partial`: as que já tinham aceitado subpassos,
`active_set_unsettled`: subpassos aceitos com o conjunto alagado/seco ainda mudando).

Outros solos, geometrias e chuvas (sem estas opções os resultados são os de antes):

| opção | significado |
|---|---|
| `c=`, `phi=`, `gamma=`, `E=`, `nu=`, `ks=`, `poro=` | resistência (kPa, graus), peso (kN/m³), elasticidade (kPa), condutividade (m/s), porosidade |
| `chi=1` \| `chi=<valor>` \| `chi=sw` \| `chi=se` | Bishop: constante só na sucção (tan φ^b = χ tan φ'; χ = 1 com p_w > p_a), S_w ou S_e |
| `swcc=fx:a,n,m` \| `swcc=vg:alpha,n,swr` | retenção de Fredlund–Xing (a em kPa) ou van Genuchten (alpha em 1/kPa) |
| `krw=pot:p` \| `krw=mualem` | permeabilidade relativa da água |
| `camadas=<arquivo>` | uma linha por camada (material id da malha): `id chave=valor ...` com as chaves acima |
| `mesh=<.msh>`, `na=x:y,x:y,...`, `lateral=0`, `basena=1` | malha qualquer (`scripts/malha_poligono.py`, camadas paralelas à superfície), NA poligonal, laterais fechadas, carga do NA na base |
| `chuva=<arquivo>` | hietograma `t0,t1,I` (h, h, mm/h) |
| `padrao=adiantado\|atrasado` | hietogramas lineares de média `I` em `train` (Tang et al. 2018) |
| `evap=` (mm/dia), `psimin=` (kPa) | evaporação quando não chove, limitada pela sucção `psimin` |
| `regime=<q mm/h>` | estado inicial: infiltração permanente q (pseudotransiente até o regime) |
| `alagado=1`, `dtmax=` | superfície alagada (p_w = 0) enquanto chove (I → ∞); passos crescentes até `dtmax` h |
| `dtchuva=` | passo durante a chuva (h, 0,25) |

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
| f10 | 180 / 360 / 900 | 0,85 / 0,81 / 0,81 | 1,00 / 1,00 / 0,99 | 0,99 / 1,00 / 0,99 | **0,86** / **0,84**\* / **0,84**\* |
| f50 | 3,6 / 9 / 18 | 1,98 / 1,79 / 1,49 | 1,96 / 1,74 / 1,32 | 2,01 / 1,86 / 1,62 | 2,01 / 1,85 / 1,61 |
| f50 | 36 / 54 / 90 / 180 | 1,05 / 0,89 / 0,82 / 0,82 | 1,01 / 1,00 / 1,00 / 1,00 | 1,21 / 1,03 / 1,00 / 0,90 | 1,18 / **0,97 / 0,85 / 0,84** |
| f100 | 0,9 / 3,6 / 9 | 2,04 / 1,88 / 1,55 | 2,01 / 1,77 / 1,36 | 2,06 / 1,95 / 1,70 | 2,06 / 1,94 / 1,69 |
| f100 | 18 / 36 | 1,45 / 1,45 | 1,28 / 1,28 | 1,52 / 1,49 | 1,50 / 1,48 |

\* refeitos depois da correção da busca linear; antes a percolação parava às 8,5 h com toda a superfície alagada
(chuva de 1 e 2,5 k_s no solo mais permeável) e o mínimo ficava superestimado (1,05 e 0,98). A coluna desacoplada
inteira foi recalculada depois das três correções do avanço no tempo (Limitações) e repetiu os demais valores.

* **Com as hipóteses do artigo** (desacoplado, ar passivo) o FS segue a sigmoide da Eq. 6 nos três solos, inclusive o
  patamar inferior (0,84–0,85 contra 0,81–0,82) e a recuperação após a chuva da Fig. 5b.
* **No modelo completo o FS mínimo não desce abaixo de ≈ 1,00.** O transiente acoplado nunca perde o equilíbrio: com
  Mohr–Coulomb associativo (ψ = φ') a plastificação dilata o solo e, na percolação acoplada, isso gera sucção
  (p_w cai) e aumenta a resistência — o talude se estabiliza no limite. Com ψ = φ'/2 o patamar é o mesmo
  (1,00–1,01); só com ψ = 0 o transiente perde o equilíbrio e o FS desce a 0,96 (contraexemplo A1, abaixo).
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

### Contraexemplos A1 e A2 (modelo completo; no A2 também ar passivo e o modelo do artigo)

`scripts/contraexemplos.py rodar --itens A1` (15 casos, executável `SlopeRainfallPV`) e `scripts/fase_ar.py`
(22 casos); figuras `contra_A1_dilatancia.png` e `fase_ar_*.png`.

**A1 — dilatância** (f50, α = 45°, acoplado com ar ativo, ψ só no transiente, FS com ψ = φ'):

| I_r (mm/h) | 18 | 36 | 54 | 90 | 180 |
|---|---|---|---|---|---|
| ψ = 0 | 1,311 | **0,962** | **0,963** | **0,960** | **0,958** |
| ψ = 13° | 1,319 | 1,005 | 1,004 | 1,005 | 1,008 |
| ψ = 26° | 1,319 | 1,003 | 1,000 | 1,002 | 1,002 |

* Sem dilatância o modelo acoplado rompe (com ψ = 0 e I_r ≥ 36 mm/h o transiente perde o equilíbrio entre 11,75
  e 15,5 h); com ψ = φ'/2 ou φ' (e I_r ≥ 36 mm/h) o FS para em 1,00–1,01. O patamar com ψ = 0 (0,96) ainda fica acima do desacoplado
  (0,84–0,85). Com 18 mm/h ψ quase não importa: a dilatância só age perto da ruptura.

**A2 — fase de ar** (f100; G1: α = 63,4°, H_s = 10 m; G2: α = 45°, H_s = 20 m):

| G1, I_r (mm/h) | 9 | 18 | 36 | 80 |
|---|---|---|---|---|
| desacoplado, p_a = 0 (artigo) | 1,38 | 1,16 | 1,10 | 1,09 |
| acoplado, ar passivo | 1,40 | 1,20 | 1,13 | 1,12 |
| acoplado, ar ativo | **1,17** | **1,03** | **1,02** | **1,01** |

* Com ar ativo o FS mínimo fica 0,08–0,22 (G1) e 0,08–0,11 (G2; no modelo do artigo o FS ainda caía às 48 h) abaixo
  do modelo do artigo, sem romper; com I_r ≥ 18 mm/h o G1 fica a 1–3 % da ruptura. No G1 infiltra menos água com ar
  ativo (10,4–11,6 m³/m contra 13,5–17,9 m³/m com ar passivo, em 24 h): o ar comprimido a 36–38 kPa eleva p_w sob a
  frente. A saída do ar (k_ra,min) só importa a partir de 1e-2.

### Contraexemplos A3–A6 (modelo do artigo: desacoplado, p_a = 0)

`scripts/contraexemplos.py` (99 casos; figuras `contra_A3_anisotropia.png`, `contra_A4_succao.png`,
`contra_A5_angulo.png`, `contra_A6_volume.png` e tabelas `contra_tabela_A*.csv`; resultados em
`resultados/contraexemplos/`). Solo f50, α = 45°, H_s = 10 m, H_w = 5 m, chuva de 24 h salvo indicação; FS mínimo
por redução de resistência.

**A3 — anisotropia e a intensidade-limite** (k_1 = 1e-4 na direção θ, k_2 = 1e-5; θ = 135° é a direção da face):

| I_r (mm/h) | 9 | 18 | 36 | 54 | 90 | 180 | 360* |
|---|---|---|---|---|---|---|---|
| isotrópico 1e-5 (B_passivo_desacoplado) | 1,85 | 1,61 | 1,18 | 0,97 | 0,85 | 0,84 | — |
| θ = 0° | 1,79 | 1,59 | 1,24 | 1,08 | 1,05 | 1,05 | 1,05 |
| θ = 45° | 1,86 | 1,63 | 1,34 | 1,19 | 0,92 | 0,87 | 0,87 |
| θ = 135° | 1,94 | 1,69 | 1,22 | 0,96 | **0,65** | **0,63** | **0,63** |

\* na primeira rodada, antes da correção da busca linear, a percolação parava às 13–24 h (toda a superfície alagada).
Todos os casos da tabela foram refeitos depois das três correções do avanço no tempo (θ = 135° com 36 mm/h passou
de 1,15 para 1,22).

* O patamar inferior da Eq. 6 (≈ 0,82) não é uma propriedade do solo: com o mesmo k_s médio vai de 0,63 a 1,05 só
  pela orientação do tensor. Com a direção mais permeável paralela à face (θ = 135°) a água infiltrada pela crista
  escoa paralela à face e o FS fica 0,2 abaixo do isotrópico; com estratificação horizontal (θ = 0°) fica 0,2 acima.
* Nas chuvas fracas (≤ 18 mm/h) a ordem se inverte (θ = 135° é o mais estável; θ = 0° fica 0,02–0,06 abaixo do
  isotrópico): as orientações inclinadas (θ = 45° e 135°, k vertical 5,5e-5 m/s) drenam a água da zona rasa; a
  estratificação horizontal, com o k vertical do isotrópico, não.
* No modelo acoplado (seção Anisotropia da permeabilidade, 36 mm/h, cálculos anteriores às correções) o pior caso
  era θ = 0 (1,08 contra 1,18 com θ = 135°); aqui, com 36 mm/h, θ = 135° e θ = 0° quase empatam (1,22 contra 1,24).
  A diferença grande (0,63 contra 1,05) só aparece no patamar das chuvas fortes (≥ 90 mm/h), que não foi calculado no
  modelo acoplado com o tensor girado.

**A4 — limite da sucção inicial** (FS inicial e FS mínimo com s_max = 75 / 150 / 300 kPa):

| H_w (m) | 2,5 | 5 | 7,5 | 10 | 15 |
|---|---|---|---|---|---|
| FS inicial | 1,80 / 1,83 / 1,83 | 2,12 / 2,20 / 2,20 | 2,35 / 2,56 / 2,56 | 2,40 / 2,91 / 2,92 | 2,40 / 3,45 / 3,63 |
| FS mínimo, 9 mm/h | 1,47 / 1,49 / 1,49 | 1,86 / 1,91 / 1,91 | 2,08 / 2,23 / 2,23 | 2,23 / 2,51 / 2,53 | 2,35 / 2,92 / 3,07 |
| FS mínimo, 36 mm/h | 0,95 / 0,98 / 0,98 | 1,17 / 1,22 / 1,22 | 1,39 / 1,49 / 1,49 | 1,49 / 1,64 / 1,64 | 1,54 / 1,83 / 1,90 |

* O teto de 75 kPa do artigo pesa sobretudo com nível d'água profundo (H_w ≥ 7,5 m): com H_w = 15 m o FS inicial
  chega a +51 % (2,40 → 3,63), mas o FS mínimo ganha menos (+23 % com 36 mm/h, 1,54 → 1,90; +31 % com 9 mm/h),
  porque, para o mesmo H_w, a queda pela chuva cresce com o limite da sucção inicial (36 mm/h: 0,86 com 75 kPa, 1,63
  com 150 e 1,73 com 300).
* O efeito é monotônico: em todos os casos um limite maior dá FS mínimo maior ou igual. (Os casos não monotônicos
  da versão anterior desta tabela eram artefato dos defeitos do avanço no tempo, ver Limitações.)

**A5 — inclinação da face** (FS inicial igual nos três solos):

| α (°) | 20 | 30 | 40 | 45 | 50 | 60 | 70 |
|---|---|---|---|---|---|---|---|
| FS inicial | 2,80 | 2,44 | 2,21 | 2,12 | 2,02 | 1,84 | 1,66 |
| f10, 80 mm/h | 1,67 | 1,44 | 1,26 | 1,20 | 1,13 | 1,00 | 0,89 |
| f50, 36 mm/h | 1,63 | 1,42 | 1,24 | 1,17 | 1,12 | 1,00 | 0,89 |
| f100, 9 mm/h | 2,43 | 2,07 | 1,79 | 1,69 | 1,59 | 1,43 | 1,30 |

* O FS mínimo cai de forma monotônica com α, numa faixa (0,74–1,13 entre 20° e 70°) comparável à separação entre os
  três pares num mesmo ângulo (0,41–0,79) e maior que a separação entre solos sob a mesma chuva na série B (0,24–0,40
  a 45°, 9–36 mm/h): a geometria não é "secundária" para o FS mínimo. A queda não é linear: em módulo, a inclinação
  local cai de 0,024–0,039 por grau (20–25°) para 0,011–0,013 (65–70°); a Tabela 7 do artigo usa uma reta por par
  (F = −0,027, −0,018 e −0,025 nos três pares). A queda causada pela chuva (FS inicial − FS mínimo) diminui com α nos
  solos permeáveis (1,13–1,16 a 20° → 0,76–0,77 a 70°) e varia pouco, sem tendência monotônica, no f100 (0,36–0,43;
  com α ≤ 30° o FS ainda caía às 48 h simuladas, e ali o mínimo é o do último instante calculado).

**A6 — mesmo volume (864 mm) em durações diferentes** (I_r = 864 mm / duração):

| duração (h) | 3 | 6 | 12 | 24 | 72 |
|---|---|---|---|---|---|
| I_r (mm/h) | 288 | 144 | 72 | 36 | 12 |
| f10 (k_s = 360 mm/h) | **1,31** | 1,38 | 1,47 | 1,57 | 1,78 |
| f50 (k_s = 36 mm/h) | 1,58 | 1,37 | 1,21 | **1,17** | 1,25 |
| f100 (k_s = 3,6 mm/h) | 1,93 | 1,84 | 1,70 | 1,48 | **1,16** |

* A duração crítica de um volume V é compatível com t* ≈ V/k_s (I_r ≈ k_s): 2,4 h no f10 (abaixo da faixa testada:
  o menor FS é o da mais curta, 3 h), 24 h no f50 (onde está o menor FS) e 240 h no f100 (acima da faixa testada: o
  menor FS é o da mais longa, 72 h). Nem "chuva curta e intensa" nem "chuva longa" é o pior caso em geral: com volume
  fixo, o pior caso é a chuva com I_r mais próximo de k_s. No f50 a infiltração vai de 22 m³/m (3 h) a 58 m³/m (72 h);
  no f10 tudo infiltra (60,5 m³/m).

## Limitações

* Patamar FS ≈ 1 do modelo acoplado com ψ = φ' (acima); a comparação com o artigo deve usar `acoplado=0 passivo=1`.
* Em três casos acoplados anteriores às correções (f100: série A com α = 26,6° e 9 mm/h; B e B_passivo com 1,8 e
  0,9 mm/h) um passo não converge nem com 5 bissecções com o FS ainda em 2,0–2,4 (mínimo 1,84–2,06): o instante é
  gravado em `t_fail` e as figuras o marcam como "transiente sem convergência" (laranja), não como ruptura.
* Três defeitos do avanço no tempo foram corrigidos depois de parte dos cálculos: a busca linear impedia a
  convergência com a correção de pressão já abaixo da tolerância (0a8bad1); depois de uma tentativa com bissecção
  que aceitava parte do passo e falhava, a nova tentativa partia do estado parcial (o estado ficava à frente do
  relógio); a infiltração dos passos subdivididos só contava o último subpasso. Juntos, os dois primeiros produziam
  saltos no transiente (A1, ψ = 0, 18 mm/h: FS mínimo 1,08 em vez de 1,31). A1–A6 e a variante desacoplada da série B
  foram recalculados com as correções; as séries A–D (variantes acoplada e com ar passivo) e o estudo de anisotropia
  (acima) são anteriores e podem conter o mesmo erro. Os resultados em `resultados/` e as figuras em `figuras/` de
  A1–A6, da variante desacoplada da série B e dos benchmarks (`resultados/benchmarks/`) são os recalculados.
* Porosidade constante nos coeficientes; sem histerese de retenção; FS por elementos finitos (SRM) e não por
  equilíbrio-limite.
* Custo: a montagem do NeoPZ (elemento multifísico de 3 espaços) domina; Pardiso (MKL) recomendado.

## Alterações na biblioteca

* `TPZMatUnsaturatedUPwPa`, `TPZUnsaturatedHydraulics` (novos, Material/Plasticity).
* `TPZSpStructMatrix`: blocos parcialmente removidos pelo filtro de equações contam só as equações ativas (antes a
  montagem abortava com um componente de um nó prescrito).
* `TPZPardisoSolver`: matrizes não simétricas sempre fatoradas (a iteração CGS com o LU anterior falhava nas
  tangentes de Newton e, com K = 1, também o fallback de fatoração).
