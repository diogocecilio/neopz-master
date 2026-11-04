# Rahardjo et al. (2007), talude de 45° com os três solos, no cluster

Série A do artigo restrita a α = 45°: solos f10,-4, f50,-5 e f100,-6, cada um com Ir = 9 mm/h, 80 mm/h e 1·ks,
chuva de 24 h, 48 h simuladas e FS de 1 em 1 h. Geometria da Fig. 2 (70 m × 30/20 m), c' = 10 kPa, φ' = φb = 26°.

| Linhas de `casos.txt` | O que rodam |
|---|---|
| 1–9 | Variante (iii), hipóteses do artigo: esqueleto rígido e p_a = 0 |
| 10–18 | Variante (i), modelo completo do `main()`: u-p_w-p_a acoplados, ar ativo |
| 19–21 | Sensibilidade de f50, 36 mm/h: chuva na projeção horizontal, Δt = 0,05 h, malha fina |

`fonte/rahardjo.cpp` inclui `fonte/main.cpp`, uma cópia do `main.cpp` do projeto, e replica o laço de chuva de `main()`.

## Como rodar

```bash
git clone --branch UnsaturatedSoil https://github.com/diogocecilio/neopz-master.git
cd neopz-master/Projects/LearningMohrCoulomb2/cluster

# carregue compilador, CMake (>= 3.14) e MKL; os nomes dos módulos variam em cada cluster
JOBS=4 ./compilar.sh                # sem MKL: USAR_MKL=OFF ./compilar.sh (skyline, mais lento)

TF_H=1 ./rodar_caso.sh 19           # teste de poucos minutos; apague resultados/ depois

# ajuste --partition e os módulos em job_slurm.sbatch; escalonador padrão: SLURM
sbatch job_slurm.sbatch
squeue -u $USER
```

### Com a VPN ligada só por alguns minutos

A compilação e os casos rodam como jobs; a VPN só precisa ficar ligada para clonar e submeter:

```bash
git clone --branch UnsaturatedSoil https://github.com/diogocecilio/neopz-master.git
cd neopz-master/Projects/LearningMohrCoulomb2/cluster
nano modulos.sh        # módulos (module load ...) e, se preciso, PARTICAO
./submeter.sh          # compila num nó e, quando terminar, roda os 21 casos; pode desligar a VPN
```

Depois, numa nova conexão: `./status.sh` (ver abaixo) e, se a compilação falhar, `tail logs/*build*.out`.
`./submeter.sh 3,7` roda só esses casos; se `bin/rahardjo` já existe, a compilação é pulada. Com
`EMAIL=<endereço> ./submeter.sh`, o SLURM avisa por e-mail quando o array termina ou falha (se o cluster enviar e-mail).

### Tempos e acompanhamento

`./status.sh` mostra, para cada linha de `casos.txt`, o estado (transiente, FS, concluído, falhou), a situação na fila do
SLURM, os passos feitos, s/passo, iterações de Newton por passo, os FS feitos, s/FS, o tempo decorrido, quanto falta no
ritmo atual e a máquina. Ele lê o que o cálculo grava a cada passo (`balanco.csv`) e a cada FS (`fs.csv`), então mostra o
progresso enquanto o job roda. Uma conexão curta basta; conectado, `watch -n 60 ./status.sh` atualiza sozinho.

```bash
./status.sh                   # tabela
./status.sh --slurm           # e o sacct do último array: estado, duração, memória máxima, nó
./status.sh --csv > tempos_fermi.csv
```

Para comparar com a sua máquina, rode lá os mesmos casos com este pacote (`./compilar.sh`, depois `./rodar_caso.sh <n>`;
`TF_H=6 ./rodar_caso.sh <n>` serve para um teste mais curto), traga os resultados do cluster e compare:

```bash
rsync -av <usuario>@<cluster>:neopz-master/Projects/LearningMohrCoulomb2/cluster/resultados/ resultados_cluster/
./status.sh --comparar resultados resultados_cluster
```

A comparação usa só os passos e os instantes de FS que as duas pastas têm, e mostra s/passo, s/FS, a razão B/A e a maior
diferença de `FS_SRM` entre as duas, que deve ser pequena se o cálculo for o mesmo.

Sem SLURM, rode `./rodar_caso.sh <n>` para n = 1 a 21, com o escalonador do cluster ou dentro de `tmux`. Não rode o
cálculo no nó de login.

## Custo

Cada tarefa usa 1 CPU, menos de 1 GB de memória e um único processo. Na estação de trabalho, com carga, foram medidos cerca
de 1,3 s por passo e 80 s por FS na variante (iii). O custo do modelo completo (i) não foi medido até o fim. A estimativa é
de cerca de 2 h por tarefa; o limite de `--time=06:00:00` é folgado.

As tarefas podem ser retomadas: submetida de novo, cada uma pula o transiente concluído e os instantes já gravados em
`fs.csv`. Para rodar só as que faltam: `./submeter.sh 3,7`. O tempo de cada passo fica em `balanco.csv` (`seg_passo`), o
de cada FS em `fs.csv` (`seg`) e o de cada fase em `tempos.csv`; `./status.sh` resume os três.

## Resultados

Ficam em `resultados/<variante>/<caso>/`:

| Arquivo | Conteúdo |
|---|---|
| `fs.csv` | FS por instante, na ordem de prioridade (t = 24, 0, 12, 48, ...); ordene por `t_h` para plotar |
| `balanco.csv` | Balanço de massa de água por passo: chuva, infiltração, laterais, armazenamento, nós alagados |
| `state_<passo>.bin` | `neq` (int64) seguido de X (double) em cada instante; passo = t/Δt |
| `nodes.csv`, `tri.csv` | Malha e equações de p_w, p_a e u em cada nó |
| `mec_pontos.scal_vec.<passo>.vtk` | Último estado convergido da busca de SRM em t = 24 h |
| `exit.txt`, `run.err`, `fs.err` | Situação do transiente e mensagens de erro |
| `tempos.csv` | Uma linha por fase executada (transiente, fs): início, fim, segundos, código de saída, máquina, job, CPU |
| `plano.txt` | Número de passos e de instantes de FS previstos para o caso |

Colunas de `fs.csv`:

| Coluna | Significado |
|---|---|
| `FS_SRM` | FS com a busca corrigida (ver abaixo); é o valor a comparar com o artigo |
| `FS_SRM_1a_parada` | Onde a busca corrigida teria parado sem a nova tentativa |
| `FS_projeto` | `FactorSRM` do projeto, sem alteração e com a lista original de F, em t = 0, 6, ..., 48 h |
| `FS_GIM`, `FS_GIM_1a_parada` | Aumento da gravidade com a busca corrigida, em t = 0, 12, ..., 48 h |
| `seg` | Tempo gasto no instante |

Valores `-1` indicam que a grandeza não foi calculada naquele instante.

Para trazer de volta: `rsync -av <cluster>:<pasta>/resultados/ <destino>/`

## Por que o FS não usa `FactorSRM` diretamente

`FactorSRM`, `FactorGIM` e `FactorLoad` em `main.cpp` tratam a primeira falha de Newton como início da bisseção e nunca
tentam de novo o valor que falhou. Quando a falha vem de um salto grande, e não do colapso, o FS fica grudado logo abaixo de
um valor da lista. Exemplo com f50, 36 mm/h, em t = 12 h, sempre sobre o mesmo campo de pressão:

| Rotina | FS |
|---|---|
| `FactorSRM` do projeto: F = 1,2 falha em 25 iterações e a busca para em 1,19375 | 1,194 |
| Mesma busca tentando F = 1,2 de novo a partir de 1,19375: converge em 5 iterações | 1,5625 |
| `SlopeRainfall` antigo, cujo `CollapseSearch` tinha essa nova tentativa | 1,565 |
| Bishop simplificado, método do artigo | 1,542 |

`fonte/rahardjo.cpp` traz `FactorSRMRetry` e `FactorGIMRetry`, com a nova tentativa e uma lista de F de 0,1 em 0,1. O
`main.cpp` do projeto não foi alterado.
