# Cam-Clay modificado no NeoPZ

Porte para C++ do modelo de `python/modified_cam_clay.py` (return mapping implícito do Cam-Clay
modificado, de Souza Neto, Perić & Owen 2008, Seção 10.1, com as extensões do RS2/FLAC3D).

## Classes

`Material/Plasticity/TPZModifiedCamClay.{h,cpp}`

* `TPZYCModifiedCamClay` – superfície de escoamento `f(p, q, a)` (eqs. 10.1-10.2) e endurecimento
  `p_c(α) = p_c0 exp(v0 α/(λ-κ))`, `a(α) = (p_c + p_t)/(1 + β)`.
* `TPZModifiedCamClay` – elasticidade + return mapping + tangente consistente 6x6. Deriva de
  `TPZPlasticBase` e pode ser usado diretamente em `TPZMatElastoPlastic<TPZModifiedCamClay, TPZElastoPlasticMem>`
  (3D) e `TPZMatElastoPlastic2D<TPZModifiedCamClay, TPZElastoPlasticMem>` (deformação plana); as duas
  instâncias já estão compiladas na biblioteca (`BUILD_PLASTICITY_MATERIALS=ON`).

| Python (`modified_cam_clay.py`)          | C++ (`TPZModifiedCamClay`)                                   |
|------------------------------------------|--------------------------------------------------------------|
| `mcc_parameters(M, lam, kap, N, v0, pc0, p0, pt, beta, elasticity, shear, G, nu, sigma0)` | `SetUp(M, lambda, kappa, N, v0, pc0, p0, pt, beta, elasticity, shear, G, nu)` + `SetInitialStress(sigma0)` (`v0 <= 0` equivale a `v0=None`) |
| `elasticity='linear' / 'pressure_dependent'` | `ELinear` / `EPressureDependent`                          |
| `shear='constant_G' / 'constant_nu' / 'hypo_nu'` | `EConstantG` / `EConstantNu` / `EHypoNu`               |
| `P['K0'] = ...` (sobrescrever K0)         | `SetLinearBulkModulus(K0)`                                   |
| `pressure`, `shear_modulus`, `hardening`, `yield_function`, `pc_of_alpha` | `Pressure`, `ShearModulus`, `YC().Hardening`, `YC().YieldValue`, `Pc` |
| `return_mapping(P, eps, epsp_n, alpha_n, sig_n=, eps_n=)` | `ReturnMapping(eps, epsp_n, alpha_n, res, &sig_n, &eps_n)` (`TResult`) |
| `ReturnMappingError`                      | `TPZModifiedCamClay::ReturnMappingError` (exceção)           |
| `invariants`                              | `TPZModifiedCamClay::Invariants`                             |

Convenções (as mesmas do Python e de `TPZMatElastoPlastic` neste ramo): tração positiva,
Voigt `{xx, xy, xz, yy, yz, zz}` (ordem de `TPZTensor`), distorções de engenharia nas deformações,
`α = -ε_v^p` guardado em `TPZPlasticState::m_hardening`, `ε^e` medida a partir do estado inicial (σ = σ0
para ε = 0).

Com `TPZMatElastoPlastic`:

* `SetPlasticityModel` cria a memória padrão com `m_sigma = σ0` (a tensão inicial faz parte do modelo; as
  cargas de contorno devem estar em equilíbrio com ela, p.ex. tração `σ0·n` nas faces livres);
* no modo `EHypoNu` a tensão do passo anterior σ_n é a `m_sigma` da memória, que o material passa como
  valor de entrada de `sigma` em `ApplyStrainComputeSigma` (um tensor nulo é lido como σ0);
* σ0 variável no espaço (p.ex. geostática) pode ser imposta sobrescrevendo
  `UpdateMaterialCoeficients(x, plasticity)` e chamando `plasticity.SetInitialStress(σ0(x))`.

## Verificação (`PlasticityTestsCamClay`)

```
cmake -DBUILD_PLASTICITY_MATERIALS=ON <neopz>
make PlasticityTestsCamClay
cd Projects2/PlasticityTestsCamClay && ./PlasticityTestsCamClay
python3 <neopz>/Projects2/PlasticityTestsCamClay/compara_python.py   # opcional (numpy, matplotlib)
```

O executável imprime cada verificação com `OK`/`FALHOU` (código de saída 1 se alguma falhar):

1. ponto material: estado OCR = 5 do RS2 x Python; tangente consistente x diferenças finitas em 200 estados
   aleatórios (p_t = 20, β = 0.6, K = -v0 p/κ, σ0 anisotrópica) e no modo hipoelástico;
2. ensaios triaxiais drenados das Figs. 8.5-8.8 do RS2 (e OCR = 5 com K dependente de p) x Python e x
   solução analítica em forma fechada, e convergência de 1ª ordem com o passo;
3. benchmark da Itasca (FLAC3D) com um hexaedro e `TPZMatElastoPlastic`: drenado e não drenado, R = 1.6 e 8,
   x elementos finitos do Python (`poro_camclay_fem.py`) e x FLAC3D. O caso não drenado usa o limite de
   permeabilidade nula do sistema u-p de Biot (pressão de poros condensada no ponto de integração,
   `p = -α M ε_v`, classe `TMatCamClayUndrained` no `main.cpp`);
4. elasticidade hipoelástica (parâmetros do benchmark 1.15.2 do Abaqus) em FE x ponto material.

`compara_python.py` recalcula os mesmos casos com o código Python de `python/` e compara curva a curva com os
CSV gravados pelo executável.
