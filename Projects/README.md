# Verificação do modelo não saturado contra a literatura

Branch `Verificacao`. Cada projeto reproduz um exemplo publicado, para verificar, validar e comparar o código: o
material `TPZMatUnsaturatedUPwPa` (u–p_w–p_a, Bishop, Mohr-Coulomb associado), a análise `TPZUnsatAnalysis`, o fator
de segurança por redução de resistência (SRM), os campos aleatórios de Karhunen-Loève e o Monte Carlo.

| projeto | referência | o que verifica |
|---|---|---|
| `V1_LiakopoulosTrifasico` | Laloui, Klubertanz & Vulliet (2003), IJNAMG 27:183-206 | coluna de Liakopoulos com ar ativo (u–p_w–p_a) e passivo (u–p_w): p_w, p_a, deslocamento e vazão contra o experimento e o benchmark ALERT |
| `V2_TaludesCho` | Liao et al. (2025), Comput. Geotech. 179:107022; Cho (2010), JGGE 136(7) | FS por SRM de dois taludes secos clássicos e P_f por Monte Carlo com campo lognormal de c_u |
| `V3_SuperficieLivre` | Sharma, Fujisawa & Murakami (2021), FEAD 197:103632 | percolação permanente com superfície livre e face de percolação (barragens retangular e de face inclinada) |
| `V4_TaludeLiu` | Liu et al. (2024), Comput. Geotech. 165:105938 | talude não saturado sob chuva permanente: poropressões, FS e Monte Carlo com campos de c', φ' e k_s |
| `V5_EncostaDuasCamadas` | Abbasov et al. (2024), Comput. Geotech. 167:106102 | encosta em duas camadas sob chuva (benchmark de Lu et al. 2012): umidade e fator de segurança local |
| `V6_MediaLocalKL` | Yuan et al. (2021), Comput. Geotech. 137:104247 | o gerador de Karhunen-Loève: a variância da média do campo ao longo de um arco |

## Resultados (resumo; detalhes no README de cada projeto)

| projeto | principal comparação | artigo | nosso | situação |
|---|---|---|---|---|
| V1 | ar passivo (regressão do LearningLiakopoulos / Callari e Abati 2009) | tabela validada | igual em todos os dígitos | ok |
| V1 | mínimo de p_a em 30 min, lei do ar do benchmark ALERT (k_ra ≥ 1e-4) | −6,95 a −5,28 kPa (5 equipes) | −5,49 kPa em 0,80 m, perfil dentro da envoltória | ok |
| V1 | o mesmo, com os dados do ar impressos (k_ra ≥ 1%) | −6,95 a −5,28 kPa | −0,58 kPa (10× menor) | não reproduz; teste `[!shouldfail]` |
| V2 | FS do talude não drenado de Cho (SRM) | 1,369 | 1,344 (−1,8%) | ok |
| V2 | FS do talude c-φ (SRM) | 1,205 (Bishop); 1,284 (MEF) | 1,219 (+1,2%); o 1,284 não é reproduzido | ok / diferença |
| V2 | P_f do talude não drenado, campo gaussiano | 10,12% (10.000 simulações) | 9,38%, IC 95% [7,5; 11,6]% (800 simulações) | ok |
| V3 | vazão da barragem retangular contra Dupuit-Charny (exato) | 1,7168e-2 m³/s/m | +0,31% (efeito da franja capilar; 0,005% com a identidade corrigida) | ok |
| V3 | ponto de saída, barragens retangular e inclinada | 12,9 m; 3,4 m | 12,5 a 13,0 m; 3,375 a 3,5 m (com a = 2/kPa) | ok |
| V4 | poropressões permanentes e nível d'água nas seções X-X e Y-Y | Fig. 3 | nível d'água ±0,02 m; sucção na crista 10% menor (permeabilidade FXH aproximada) | ok, com ressalva |
| V4 | FS determinístico (SRM contra Bishop) | 1,525 | 1,588 (+4,1%) | ok |
| V5 | mapas de umidade, sucção e LFS (Figs. 4 a 6) | mapas | RMS de θ 0,004 a 0,008; perfis do LFS a 0,02 a 0,05 | ok |
| V5 | área com LFS < 1 em 15 h | 14,9 m² | 22,8 m² (+53%) | diferença aberta |
| V6 | desvio da média do campo KL num arco contra a Eq. 14 | Eq. 14 | −0,23% a −0,30% (truncamento de 90%), Monte Carlo dentro de 1,9 erros-padrão | ok |
| V6 | Tabelas 5 e 6 (σ da média no arco) | valores impressos | 0 a 19% menores (os arcos do artigo não são dados) | não reproduz |

## Rodar no cluster

O pacote está em `Cluster/` (ver `Cluster/README.md`): compilação de tudo no nó de login, os testes num nó, os
casos completos e o Monte Carlo de V2 e V4. Exemplos, de dentro de `Projects/Cluster` no fermi:

```bash
source modulos.sh && ./compilar.sh
./submeter_testes.sh
./submeter_caso.sh V3 retangular inclinada
mc/submeter_mc.sh V2 4
./status.sh
```

Cada projeto tem também uma pasta `cluster/` com scripts próprios, mais simples (compilar e um sbatch). Prefira o
pacote `Cluster/`, que compila uma vez para todos.

## Estrutura

- `Comum/`: código compartilhado, só cabeçalhos.
  - `UnsatModel.h`, `UnsatAnalysis.h`, `UnsatProblem.h`, `UnsatFS.h`, `UnsatReference.h` e `UnsatTaylor.h`: o
    modelo, as malhas estruturadas, a análise e os fatores de segurança.
  - `RahardjoDriver.h`: parâmetros do solo, leitura de malha Gmsh, malhas computacionais, `FactorSRM`, `FactorGIM`,
    `FactorLoad` e VTK. É o antigo driver do talude de Rahardjo; inclua-o num único `.cpp` por executável.
  - `KLSeparavel.h`: campos de Karhunen-Loève com covariância exponencial separável.
  - `malhas/`: as malhas do talude de Rahardjo.
- `Cluster/`: scripts para o fermi (SLURM).
  - `modulos.sh`: os módulos.
  - `mc/`: a fila de Monte Carlo e o pós-processamento herdados do projeto `SlopeRandomField`.
- `V1_...` a `V6_...`: um projeto por verificação. Cada um tem:
  - `problema.h`: parâmetros, malha, condições de contorno e as funções que rodam os casos;
  - `main.cpp`: o executável do caso completo;
  - `testes.cpp`: os testes Catch2 registrados no CTest;
  - `README.md`: o caso, os dados do artigo, os resultados e a comparação;
  - os dados de referência digitalizados.

## Compilar e testar

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON
cmake --build build -j 4
cd build/Projects && ctest
```

Os testes demorados têm a marca `[.lento]` e ficam fora do `ctest` padrão. Para rodá-los:

```bash
./V1_LiakopoulosTrifasico/V1_testes "[.lento]"
```

Os projetos de antes (Learning*, Unsat*, SlopeRandomField) continuam no branch `UnsaturatedSoil`.
