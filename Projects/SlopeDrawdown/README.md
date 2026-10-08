# SlopeDrawdown

Rebaixamento rápido do reservatório em frente ao talude de `SlopeMohrCoulomb` (solo saturado, nível
d'água inicial no topo) e fator de segurança (FS) em cada instante, por **aumento de gravidade** e
por **redução de resistência (SRM)** — o mesmo procedimento de `SlopeMohrCoulomb`.

Cenário de referência: Ceron, Cecílio, Linn & Maghous, *Stability analysis of slope subjected to
seepage forces considering spatial variability of soil properties*, IJNAMG 49(11), 2025
(doi:10.1002/nag.3993): c = 10 kPa, φ = 30°, γ = 20 kN/m³, rebaixamento total (h_w = H), forças de
percolação = −∇p como forças de volume, resistência em tensões efetivas.

## Modelo

1. **Análise u-p (Biot) poro-elastoplástica** — `TPZMatPoroElastoPlasticUP` com o Mohr–Coulomb do
   artigo (`TPZPlasticStepVoigt<TPZYCMohrCoulombPV2>`) e o driver `TPZPoroElastoPlasticUPAnalysis`;
   Taylor–Hood (u P2, p P1) em `TriGMesh(2)`.
   * Contornos mecânicos: base fixa, laterais em rolete. Hidráulicos: base e laterais impermeáveis;
     superfície do terreno (pé −3, topo −4, face −6) com p = γw (H − y)⁺ e carga d'água t = −γw (H − y)⁺ n
     nas partes submersas (elementos gêmeos −13 e −16, pois um id de contorno tem uma só condição).
     Acima d'água a superfície é drenada (p = 0): a linha freática permanece no topo, como no
     rebaixamento rápido de Ceron et al.
   * O nível H = H0 − λ(H0 − H1) é lido do fator de carga λ do passo (a bissecção do driver o interpola).
   * Passos: (a) gravidade com o reservatório no topo em um passo drenado (regime permanente);
     (b) rebaixamento instantâneo: passo não drenado (Δt = 0, λ: 0 → 1); (c) adensamento
     (percolação transiente), T = c_v t / H² = 10⁻³ … 10 (4 passos por década) e regime permanente.
2. **Estabilidade** em cada estado (`dry`, `crest`, `T0`, `T0.1`, `T1`, `T10`, `steady`):
   `SlopeAnalysis.h` sem alteração, malha elastoplástica de `SlopeMohrCoulomb` com a força de volume
   substituída por b = λ(γ_sat g − ∇p⁺), p⁺ = max(p, 0) (sucção desprezada), p congelado do estado.
   * No aumento de gravidade λ multiplica γ_sat **e** ∇p (pesos do solo e da água), o que equivale a
     c/λ; no SRM λ = 1 e c, tan φ são reduzidos.
   * Com p = p_água no contorno submerso, a tração efetiva é nula: o problema efetivo com b é
     equivalente ao problema em tensões totais com a carga d'água.

Funções pedidas (u-p): `CreateAtomicMesh` (o `CreatePressureMesh` original generalizado: `nstate = 1`
pressão, `2` deslocamento, com **todos** os ids de contorno nas duas malhas) e `CreateCompMesh`
(material, contornos e espaço multifísico com memória).

Parâmetros: E = 20000 kPa, ν = 0,3 (drenado do esqueleto; a resposta não drenada vem do
acoplamento), c = 10 kPa, φ = ψ = 30°, γ_sat = 20 kN/m³, γw = 10 kN/m³, Biot α = 1, 1/M = 0
(grãos e água incompressíveis), k = 10⁻² m/dia (só muda a escala de tempo: os estados em T não
dependem de k).

## Compilar e executar

```
ninja SlopeDrawdown
./SlopeDrawdown            # 2 ciclos de refinamento da zona plástica
./SlopeDrawdown nref=3
./SlopeDrawdown L=35       # rebaixamento parcial (nível final y = 35)
```

Saídas: `drawdown_up.scal_vec.N.vtk` (pressão e deslocamento da análise u-p, N = estado) e
`drawdown_<estado>.scal_vec.0.vtk` (variáveis plásticas no colapso do SRM).

## Resultados

FS com a malha inicial (`nref=0`; o refinamento da zona plástica reduz o FS, ver `SlopeMohrCoulomb`)
e Bishop simplificado com o mesmo campo p⁺ (busca de círculos, script independente):

| estado | T = c_v t / H² | FS gravidade | FS SRM | Bishop (SRM) |
|---|---|---|---|---|
| seco | — | 3,041 | 1,399 | 1,207 |
| reservatório no topo | — | 6,060 | 1,924 | 1,604 |
| logo após o rebaixamento (não drenado) | 0 | 1,544 | 1,137 | 1,038 |
| adensamento | 0,1 | 1,646 | 1,164 | |
| adensamento | 1 | 1,358 | 1,141 | 0,933 |
| adensamento | 10 | 1,288 | 1,125 | |
| regime permanente (campo de Ceron et al.) | ∞ | 1,235 | 1,106 | 0,885 |

## Limitações

* Modelo saturado: a face acima d'água é tratada com p = 0 (face de percolação) e a sucção é
  desprezada no FS (p⁺); o rebaixamento da linha freática no maciço (fluxo não confinado) não é
  modelado.
* `TPZPlasticStepVoigt` é formulado em deformação total: o estado inicial vem do passo de gravidade
  (σ' inicial da memória é ignorado).
* Mohr–Coulomb associativo: no passo não drenado a dilatância plástica geraria sucção; o acréscimo
  de resistência correspondente é desprezado no FS (p⁺, análise drenada com p congelado).
