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

O caso `resultados/B/f50_a45_H10_Hw5_I36` traz os arquivos completos e a animação.

RESULTADOS
