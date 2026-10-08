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
   Taylor–Hood (u P2, p P1) em `TriGMesh(ref)` (`ref=1`, a malha inicial de `SlopeMohrCoulomb`).
   * Contornos mecânicos: base fixa, laterais em rolete. Hidráulicos: base e laterais impermeáveis;
     superfície do terreno (pé −3, topo −4, face −6) com p = γw (H − y)⁺ e carga d'água t = −γw (H − y)⁺ n
     nas partes submersas (elementos gêmeos −13 e −16, pois um id de contorno tem uma só condição).
     Acima d'água a superfície é drenada (p = 0): a linha freática permanece no topo, como no
     rebaixamento rápido de Ceron et al.
   * O nível H = H0 − λ(H0 − H1) é lido do fator de carga λ do passo (a bissecção do driver o interpola).
   * Passos: (a) gravidade com o reservatório no topo em um passo drenado (regime permanente);
     (b) rebaixamento instantâneo: passo não drenado (Δt = 0, λ: 0 → 1); (c) adensamento
     (percolação transiente), T = c_v t / H² = 10⁻³ … 10 (4 passos por década) e regime permanente.
2. **Estabilidade** em cada estado (`dry`, `crest`, `T0`, `T0.1`, `T1`, `T10`, `steady` e a referência
   `drained`: linha freática no nível final, hidrostática, sem recarga pelo topo):
   `SlopeAnalysis.h` sem alteração, malha elastoplástica de `SlopeMohrCoulomb` partindo do mesmo nível
   `TriGMesh(ref)` da malha u-p (∇p constante em cada elemento, inclusive após o refinamento), com a força
   de volume substituída por b = λ(γ_sat g − ∇p⁺), p⁺ = max(p, 0) (sucção desprezada), p congelado.
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
./SlopeDrawdown ref=2      # malhas u-p e de estabilidade mais finas (bem mais lento)
./SlopeDrawdown L=35       # rebaixamento parcial (nível final y = 35; L em [30, 40])
```

Saídas: `drawdown_up.scal_vec.N.vtk` (pressão e deslocamento da análise u-p, N = estado) e
`drawdown_<estado>.scal_vec.0.vtk` (variáveis plásticas no colapso do SRM).

## Resultados

`./SlopeDrawdown nref=3` (malhas `TriGMesh(1)`, 32 min). FS no ciclo 3 e Bishop simplificado com o mesmo
campo p⁺ (busca multi-início, `docs/slope_mohr_coulomb/scripts/bishop_pw.py`):

| estado | T = c_v t / H² | FS gravidade | FS SRM | Bishop gravidade | Bishop SRM |
|---|---|---|---|---|---|
| seco | — | 1,906 | 1,229 | 1,806 | 1,206 |
| reservatório no topo | — | 3,799 | 1,637 | 3,611 | 1,597 |
| logo após o rebaixamento (não drenado) | 0 | 1,000 | 1,000 | 0,943 | 0,974 |
| adensamento | 0,1 | 0,910 | 0,959 | 0,854 | 0,931 |
| adensamento | 1 | 0,793 | 0,895 | 0,737 | 0,860 |
| adensamento | 10 | 0,770 | 0,881 | 0,718 | 0,846 |
| percolação permanente (campo de Ceron et al.) | ∞ | 0,781 | 0,886 | 0,723 | 0,850 |
| drenado (linha freática no pé) | — | 1,906 | 1,230 | 1,806 | 1,206 |

* O rebaixamento rápido leva o talude ao equilíbrio-limite (FS = 1,00) e, com a linha freática mantida
  no topo, a percolação o torna instável (FS ≈ 0,88 no regime permanente; Γ = 20·FS_grav ≈ 15,6 contra
  ≈ 36 do talude seco).
* Elementos finitos × Bishop: 2–4 % (SRM) e 5–8 % (gravidade) acima no ciclo 3, o mesmo padrão do
  talude seco, que converge para Bishop com mais refinamento.
* Com os campos de `ref=2` Bishop dá FS 4–8 % maior (1,038 logo após o rebaixamento, 0,884 no regime
  permanente): a poropressão junto ao pé pede malha u-p mais fina.
* Verificação da percolação: solução P1 independente (`docs/slope_mohr_coulomb/scripts/laplace_check.py`)
  igual à do NeoPZ até 5·10⁻⁴ kPa, descontado o termo de armazenamento do último passo.
* Relatório completo: `docs/slope_mohr_coulomb/relatorio.tex`.

## Limitações

* Modelo saturado: a face acima d'água é tratada com p = 0 (face de percolação) e a sucção é
  desprezada no FS (p⁺); o rebaixamento da linha freática no maciço (fluxo não confinado) não é
  modelado.
* `TPZPlasticStepVoigt` é formulado em deformação total: o estado inicial vem do passo de gravidade
  (σ' inicial da memória é ignorado).
* Mohr–Coulomb associativo (o retorno fechado exige ψ = φ): no passo não drenado a dilatância
  plástica reduz p (≈ −7 kPa em profundidade, contra ≈ −28 kPa da variação causada pelo
  rebaixamento). Na malha inicial, com o modelo antigo e ψ = 5°, o FS logo após o rebaixamento é
  1,106 (SRM) e 1,370 (gravidade), contra 1,137 e 1,540 com ψ = φ: esse FS é otimista em alguns por
  cento. O regime permanente não depende de ψ.
* A análise u-p não é refinada: com `ref=1` ela não rompe (FS do ciclo 0 > 1); com `ref=2` uma cunha
  junto à face se desloca sem limite no último passo (FS < 1). A poropressão continua sendo a da
  percolação permanente (diferença ≤ 3 kPa), mas esses deslocamentos não têm significado físico.
* O campo de poropressão depende da malha u-p: Bishop com o campo de `ref=2` dá FS cerca de 4–6 %
  maior que com o de `ref=1`.
* Com o Mohr–Coulomb de deformação total, uma tensão efetiva inicial dada em `InitializeMemory` é
  convertida na deformação própria εᵖ = ε − Dᵉ⁻¹σ'₀ (correção feita no material u-p).
