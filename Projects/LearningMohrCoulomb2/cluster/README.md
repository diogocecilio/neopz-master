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
No fermi os nós de computação não têm `cmake` (o de login tem), então `modulos.sh` define `COMPILAR=login`: `submeter.sh`
compila no nó de login em segundo plano (`nohup`, ~10 min com 4 núcleos, log em `logs/compilar_login.log`) e, ao terminar,
submete os casos sozinho; pode desconectar logo depois. O executável compilado com o gcc 14 do login roda nos nós, que têm
gcc 15, e o job de casos confere com `ldd`, antes de começar, que o nó acha todas as bibliotecas. Em outro cluster,
`COMPILAR=job ./submeter.sh` compila num job (`job_compilar.sbatch`), do qual os casos dependem.

Os casos vão num único job (`job_casos.sbatch`), que conta como um job na cota do SLURM: no fermi, um array de 21 tarefas
excede a cota de jobs submetidos por usuário (`QOSMaxSubmitJobPerUserLimit`). Esse job pede um nó com `NCPU` CPUs (padrão:
uma por caso, limitado aos núcleos do nó) e roda até `NCPU` casos ao mesmo tempo, um processo de 1 CPU cada; quando um termina,
o próximo da lista começa. `NCPU=8 ./submeter.sh` roda 8 por vez, `TEMPO=2-00:00:00` muda o limite de tempo (padrão 1 dia;
os casos são retomáveis) e `MODO=array` volta ao array de `job_slurm.sbatch`.

`./submeter.sh 3,7` roda só esses casos; se `bin/rahardjo` já existe, a compilação é pulada, e se a compilação já está na
fila, os casos esperam por ela. Com outro job de casos na fila, `submeter.sh` não submete (as duas cópias gravariam nas
mesmas pastas), a menos que `FORCAR=1`. Com `EMAIL=<endereço> ./submeter.sh`, o SLURM avisa por e-mail quando os casos
terminam ou falham (se o cluster enviar e-mail). A saída de cada caso fica em `logs/casos/<n>.out`.

### Tempos e acompanhamento

`./status.sh` mostra, para cada linha de `casos.txt`, o estado (transiente, FS, concluído, falhou), a situação na fila do
SLURM, os passos feitos, s/passo, iterações de Newton por passo, os FS feitos, s/FS, o tempo decorrido, quanto falta no
ritmo atual e a máquina. Ele lê o que o cálculo grava a cada passo (`balanco.csv`) e a cada FS (`fs.csv`), então mostra o
progresso enquanto o job roda. Uma conexão curta basta; conectado, `watch -n 60 ./status.sh` atualiza sozinho.

```bash
./status.sh                   # tabela
./status.sh --slurm           # e o sacct do último job de casos (no fermi o sacct não responde)
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

Os tempos de uma rodada com os 21 casos juntos são os de um nó cheio. Se o nó tiver hyper-threading (`lscpu`: `Thread(s) per
core: 2`), 21 processos em 22 CPUs ficam dois por núcleo físico, cada um mais lento do que sozinho. Para medir um caso como
na sua máquina, rode-o à parte, noutra pasta, com um núcleo físico só para ele:

```bash
EXCLUSIVO=1 NCPU=4 PARALELO=2 SAIDA=resultados_bench TF_H=6 FORCAR=1 ./submeter.sh 1,10
./status.sh resultados_bench
```

`EXCLUSIVO=1` reserva o nó inteiro (nada divide CPU, cache e turbo), e os dois casos, um da variante (iii) e um do modelo
completo (i), ficam cada um num núcleo físico. No fermi (Xeon Silver 4116, 12 núcleos, 2 threads por núcleo, até 3,0 GHz),
com os 21 casos juntos o clock médio ficou em 62% do máximo.

`PARALELO` é quantos casos rodam ao mesmo tempo (padrão: `NCPU`), `SAIDA` a pasta de resultados e `TF_H` as horas simuladas.
Cada caso grava em `maquina.txt` o nó e a CPU em que rodou.

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
| `maquina.txt`, `ativo` | Nó, job, início e CPU da execução; `ativo` existe enquanto o caso roda |

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
