# Solo não saturado: projetos de verificação (branch `UnsaturatedSoil`)

Modelo acoplado u–p_w–p_a com Mohr–Coulomb associado PV2 (Lira Cecílio, J. Eng. Math. 157:10, 2026):
`TPZMatUnsaturatedUPwPa<TPZPlasticStepVoigt<TPZYCMohrCoulombPV2, TPZElasticResponse>, TPZUnsaturatedMem>`.
A formulação, os resíduos, o Jacobiano e os resultados estão em
`docs/biot_nao_saturado/formulacao/formulacao.tex`. Os projetos antigos foram apagados neste branch; cada parte do
plano de verificação tem agora o seu projeto, com testes Catch2 registrados no CTest.

| projeto | parte do plano | o que verifica | tempo |
|---|---|---|---|
| `UnsatMaterialPoint` | ponto material | leis hidráulicas, Bishop com χ = S_w, PV2 contra o artigo e o ponto mais próximo, tangente (Taylor), redução de resistência | < 1 s |
| `UnsatWaterPhase` | teste 2 | água com ar passivo e esqueleto rígido: Jacobiano, hidrostático, Gardner, balanço de massa, infiltração de Philip (Callari e Abati 2009) | 3 min |
| `UnsatAirPhase` | teste 3 | água e ar com esqueleto rígido: Jacobiano, Boyle, limite de ar passivo, balanços | 3 s |
| `UnsatCoupled` | testes 0 e 4 | Jacobiano completo com plasticidade, Terzaghi, Mandel, Liakopoulos | 3 min |
| `UnsatSRM` | teste 1a | redução de resistência no talude seco contra Bishop | 4 min |
| `UnsatGIM` | teste 1b | aumento de gravidade no talude seco contra Bishop; efeito da sucção | 6 min |

`UnsatCommon/` reúne o código comum (só cabeçalhos):

* `UnsatModel.h`: parâmetros (`UnsatParams`), material, malhas estruturadas (retângulo, talude) e espaço P2/P1/P1;
* `UnsatAnalysis.h`: Newton com filtro de Dirichlet, busca linear, resíduo e produto Jacobiano–vetor;
* `UnsatFS.h`: busca do colapso por bisseção, estado geostático, SRM e GIM com pressões congeladas;
* `UnsatProblem.h`: malha + análise de um teste;
* `UnsatReference.h`: soluções de referência (Terzaghi, Gardner, Mandel);
* `UnsatTaylor.h`: teste de Taylor do Jacobiano.

## Decisões aplicadas no código

1. χ = S_w é o padrão do material (`TPZMatUnsaturatedUPwPa.h`).
2. O armazenamento do ar usa n ΔS_a + n S_a^n Δp_a / P (antes S_a^{n+1}): a massa de ar se conserva exatamente.
3. O Mohr–Coulomb é o PV2 do artigo com ψ = φ.

## Como rodar

Configurar o NeoPZ com `-DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON` (o Catch2 v2.13.10 é baixado se não for
encontrado). Na pasta de compilação:

```sh
cd Projects
ctest                                   # os seis projetos, testes rápidos
./UnsatCoupled/UnsatCoupled "Teste 4b*" # um teste pelo nome
./UnsatSRM/UnsatSRM "[.lento]"          # testes lentos (SRM e GIM com nf = 80, Philip refinado, Liakopoulos com ar ativo)
```

`UnsatCoupled` lê os dados de Liakopoulos de `UnsatCoupled/data/` pelo caminho de origem (`UNSAT_SOURCE_DIR`).

## Resultados (30/10/2026)

Todos os testes rápidos passam. Detalhes no README de cada projeto e na Seção 22 do documento.

* Jacobiano: ordem 2,000 no teste de Taylor na água, no ar e no modelo completo com 96 pontos plásticos.
* Conservação: erro do balanço de água 2·10⁻¹⁴; massa de ar 8·10⁻¹⁶ (Boyle) e 5·10⁻¹⁶ (ar livre).
* Infiltração de Philip (Callari e Abati 2009): frente a até 2,0 mm da solução semianalítica (elemento de 5 mm).
* Terzaghi e Mandel: erro máximo de p/p₀ 0,0028 e 0,0038 (centro).
* Liakopoulos (ar passivo): RMS da vazão 0,225·10⁻⁶ m/s e de p_w 0,44 kPa, contra 0,227 e 0,48 da solução de
  Callari e Abati (2009).
* Talude seco, malha nf = 80 (52 422 equações): FS_SRM = 1,2109 e FS_GIM = 1,8146, contra 1,206 e 1,806 do Bishop
  simplificado.
