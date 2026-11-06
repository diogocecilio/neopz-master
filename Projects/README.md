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
