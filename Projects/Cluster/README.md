# Projects/Cluster: os projetos de verificação no cluster (fermi)

Este pacote roda os seis projetos de verificação (V1 a V6, ver `../README.md`) no fermi (CESUP/UFRGS, SLURM): compila
tudo de uma vez no nó de login, roda os testes num nó de computação, roda os casos completos de cada projeto e o Monte
Carlo dos projetos que o têm (V2_TaludesCho e V4_TaludeLiu). Tudo volta em segundos e pode ser acompanhado de longe:
os scripts que precisam compilar antes compilam em segundo plano e submetem os jobs quando a compilação termina.

| Arquivo | O que faz |
|---|---|
| `modulos.sh` | os módulos do fermi (Lmod): o MKL do oneAPI 2025.1; os scripts o carregam sozinhos |
| `compilar.sh` | configura e compila a biblioteca e os alvos de V1 a V6 (executável e testes) em `build/`; copia os executáveis para `bin/`, com `bin/versao.txt` |
| `submeter_testes.sh`, `testes.sbatch` | os testes: o `ctest` padrão e os `[.lento]` de cada projeto, num nó de computação |
| `submeter_caso.sh`, `caso.sbatch` | os casos completos (o executável de cada projeto), em `$DATA/verificacao/<projeto>/<caso>/` |
| `mc/` | o Monte Carlo de V2 e V4: a fila (`rodar_mc.sh`, `job_mc.sbatch`, `submeter_mc.sh`), o estado (`status_mc.sh`) e o pós-processamento |
| `status.sh` | o estado de tudo: compilação, executáveis, jobs, avisos, casos completos, último teste, campanhas de Monte Carlo |
| `comum.sh` | funções comuns aos scripts (não se roda sozinho) |

Os jobs se chamam `verif-testes`, `verif-caso-<Vn>` e `verif-mc-<V2|V4>`; cada um ocupa um nó inteiro (24 CPUs, 12
núcleos físicos no fermi) e roda 12 processos de 1 CPU (com hyper-threading, dois processos por núcleo deixam cada um
~2× mais lento). As saídas dos jobs vão para `logs/<nome>_<job>.out`.

## 1. O código no fermi

```bash
ssh <usuario>@fermi.cesup.ufrgs.br
cd ~/neopz-master && git fetch origin && git checkout Verificacao && git pull origin Verificacao
cd ~/neopz-master/Projects/Cluster
```

Todos os comandos abaixo são dados nesta pasta, `Projects/Cluster`.

## 2. Módulos

```bash
source modulos.sh
module list            # deve mostrar o mkl/2025.1 (e o compilador/oneAPI que ele exige)
```

Os `submeter_*.sh` e os jobs carregam `modulos.sh` sozinhos. À mão, só antes de `./compilar.sh`. Em outro
cluster, ajuste o módulo em `modulos.sh` (`module avail`); sem MKL, `export USAR_MKL=OFF` (solver skyline, mais lento).

## 3. Compilar (no nó de login)

Os nós de computação do fermi não têm cmake: compila-se no nó de login. Em geral não é preciso chamar `compilar.sh`:
`submeter_testes.sh`, `submeter_caso.sh` e `mc/submeter_mc.sh` compilam sozinhos, em segundo plano, quando o executável
falta ou é mais velho que os fontes (depois de um `git pull`), e submetem ao fim. À mão:

```bash
source modulos.sh
./compilar.sh                  # a biblioteca e V1..V6: uns 15 a 30 min da primeira vez, alguns minutos depois
./compilar.sh V2 V4            # so alguns projetos (a biblioteca entra sempre; ver abaixo)
nohup ./compilar.sh > logs/compilar_mao.log 2>&1 &   # para poder desconectar
```

- `build/` tem a biblioteca (`libpz.so`) e, em `build/Projects/<projeto>/`, os executáveis e os testes. Os executáveis
  são copiados para `bin/` (troca atômica: uma realização que já está rodando continua com o antigo).
- `bin/versao.txt` tem uma linha por executável: `<executável> <commit>[+mod] <data>`. `+mod` quer dizer que havia
  fontes alterados ou novos fora do git. O Monte Carlo grava essa versão no `resumo.csv` de cada realização.
- Um projeto que não compila não impede os outros: fica de fora, o erro do compilador fica no log, e `compilar.sh`
  termina com `Falhou: <projeto>` e código 1. Se a biblioteca mudou nessa compilação, o executável antigo dele sai de
  `bin/` (vai para `bin/antigos/`), porque rodaria com uma `libpz.so` de outra versão. Pelo mesmo motivo, se a
  biblioteca muda ao compilar só alguns projetos, os outros que estão em `bin/` são compilados também.
- Os executáveis de `bin/` carregam a `libpz.so` de `build/`. Por isso `compilar.sh` **recusa** compilar enquanto
  houver jobs `verif-*` na fila (rodando ou esperando), processos de `bin/` ou de `build/` rodando nesta máquina, ou se
  o `squeue` não responde. `FORCAR=1 ./compilar.sh` passa por cima (os jobs em andamento ficariam com a biblioteca
  nova). Duas compilações ao mesmo tempo também são recusadas.
- Variáveis: `JOBS` (processos do make, padrão 4), `BUILD` (padrão `./build`), `USAR_MKL` (ON/OFF),
  `MKL_THREAD_MODEL` (padrão `seq`), `BAIXAR_CATCH2=0` (não baixa o Catch2 para `deps/`; o CMake o baixa ao
  configurar).

## 4. Testes (num nó de computação)

```bash
./submeter_testes.sh                    # todos: os padrao (ctest) e os [.lento] de V1..V6
./submeter_testes.sh V2 V4              # so alguns projetos
LENTOS=0 ./submeter_testes.sh           # so os padrao
V2_NMC=1000 ./submeter_testes.sh V2     # variaveis dos proprios testes (ver o README de cada projeto)
cat logs/testes_<job>/resumo.txt        # o resultado
```

- Primeiro o `ctest` em `build/Projects` (os testes padrão de cada projeto, todos ao mesmo tempo), depois os `[.lento]`
  de cada projeto (`<projeto>_testes "[.lento]"`), um processo por projeto, todos ao mesmo tempo. O tempo total é o do
  mais demorado (aqui, num núcleo: V1 ~5 min, V2 ~22 min com `V2_NMC=200`, V3 ~4 min, V4 ~10,5 min, V6 ~4 min; o V5 roda
  o caso completo na malha fina e é o mais longo). Limite do job: 12 h (`TEMPO=1-00:00:00 ./submeter_testes.sh`).
- `logs/testes_<job>/resumo.txt`: uma linha por grupo (passou ou FALHOU, o tempo, a última linha do Catch2), também na
  saída do job, `logs/verif-testes_<job>.out`. Ao lado: `ctest.log` e `<projeto>_lento.log`. O código de saída do job
  é 0 só se tudo passou.
- Os nós do fermi não têm cmake, e o `ctest` vem com ele. O job usa o `ctest` do PATH ou o do nó de login (o caminho que
  `compilar.sh` anota em `build/ctest_caminho`); se nenhum rodar no nó, roda os testes registrados no CTest direto,
  como o `ctest` os rodaria (cada executável na sua pasta de build, sem argumentos), e o `resumo.txt` diz isso.
- Os `[.lento]` rodam em `$DATA/verificacao/testes/<job>/<projeto>/` (alguns gravam arquivos na pasta corrente).
- Testes mais velhos que os fontes geram um aviso no `resumo.txt` (o job testa o que foi compilado).
- À mão: `sbatch testes.sbatch` (de dentro de `Projects/Cluster`), com `PROJETOS="V1 V3"`, `LENTOS=0`, `TIMEOUT_TESTE`.

## 5. Casos completos

Cada projeto tem um executável que roda o caso completo (a malha e os tempos do artigo), configurado por variáveis de
ambiente, e grava CSVs na pasta corrente, com `resultado.txt` por último. Uma entrada é `<caso>[:VAR=valor,...]`:

```bash
./submeter_caso.sh V1 todos
./submeter_caso.sh V3 retangular inclinada dois_blocos
./submeter_caso.sh V3 retangular_a2:CASO=retangular,AG=2 retangular_fina:CASO=retangular,NX=64,NY=96
./submeter_caso.sh V5 tabela1:VTK=1,NU_LISTA=0.2,0.25,0.3,0.35,0.4 tabela4 fig9 tabela4_ks001:CASO=tabela4,KS_INF=0.01
./submeter_caso.sh V6 tabelas b40_h10:N=100000,SEMENTE=101 fig15:N=100000,SEMENTE=110
./submeter_caso.sh V4 det                 # o deterministico do V4 (o Monte Carlo esta na secao 6)
```

- `<caso>` é o nome da pasta, `$DATA/verificacao/<projeto>/<caso>/` (`SAIDA_CASOS` troca `$DATA/verificacao`), e o
  `CASO` do programa, a menos que as trocas definam `CASO`. Os casos de cada projeto: V1 `bifasico`, `trifasico`,
  `trifasico_laloui`, `trifasico_alert` ou `todos`; V3 `retangular`, `inclinada`, `dois_blocos`; V5 `tabela1`,
  `tabela4`, `fig9`; V6 `b40_h10` ... `b40_h50`, `h50_b30` ... `h50_b50`, `fig15`, `tabelas`. As variáveis e os casos
  recomendados estão no README de cada projeto.
- Um valor pode ser uma lista com vírgulas: um pedaço sem `=` continua o valor anterior (`NU_LISTA=0.2,0.25,0.3`).
- Os casos de um pedido vão num job só (um nó), ao mesmo tempo, até 12. Limite do job: 1 dia (`TEMPO=...`).
- Na pasta: `log.txt` (a saída do programa), `caso.txt` (entrada, variáveis, versão do executável, job, início e código
  de saída), `variaveis.txt` e os CSVs do programa.
- Retomável: uma pasta com `resultado.txt` é pulada, então submeter de novo completa o que faltou. Uma pasta rodada com
  outras variáveis é recusada (dê outro nome ao caso); `REFAZER=1` roda de novo (a pasta antiga vira
  `<caso>.ant_<data>`). Uma entrada que outro job está rodando é pulada.
- Variáveis do ambiente também chegam ao programa (`T_FIM=7200 ./submeter_caso.sh V1 todos`), mas só as da entrada
  ficam em `variaveis.txt`: prefira pô-las na entrada.
- À mão: `PROJETO=V3 CASOS="retangular inclinada" sbatch caso.sbatch` (de dentro de `Projects/Cluster`).

## 6. Monte Carlo (V2 e V4)

```bash
mc/submeter_mc.sh V2 4                                       # a campanha padrao do V2 em 4 nos
CASOS="cho1 cho1:NUCLEO=exp" ULTIMA=2000 mc/submeter_mc.sh V2 2
CASOS="cho2:FS_MAX=1" mc/submeter_mc.sh V2 2                  # so P_f (busca do FS censurada em 1)
CASOS="cho1 cho1:FS_METODO=arco" mc/submeter_mc.sh V2 2       # bissecao x SRM por arc-length, mesmas realizacoes
mc/submeter_mc.sh V4 4                                       # s1, s2 e s3 com ANGULO 0 e 90
CASOS="s3 s3:RHO=-0.5 s2:ANGULO=90,RHO=-0.5" mc/submeter_mc.sh V4 2
```

O primeiro argumento é o projeto (`V2`, `V4` ou o nome completo), o segundo o número de nós (jobs), padrão 4. Os jobs
dividem a mesma fila: o caso determinístico de cada entrada (`DET=1`, o padrão) e depois as realizações `PRIMEIRA` a
`ULTIMA` de todas as entradas, intercaladas (a 1 de todas, depois a 2...), para que uma campanha cortada no meio tenha o
mesmo número de realizações em cada entrada. Rodar de novo acrescenta nós à mesma fila (nada concluído ou rodando se
repete).

| Projeto | `CASO` | Campanha padrão (`CASOS`) | Realizações | Tempo por realização (1 CPU) |
|---|---|---|---|---|
| V2_TaludesCho | `cho1`, `cho2` | `cho1 cho2` | 1 a 10 000 | ~1 min (h = 0,5 m, busca completa); ~15 s com `FS_MAX=1` |
| V4_TaludeLiu | `s1`, `s2`, `s3` | `s1 s2 s3 s1:ANGULO=90 s2:ANGULO=90 s3:ANGULO=90` | 1 a 2 000 | ~1,5 min (s1, s3); ~1 min (s2) |

- **Entradas:** `CASOS="<caso>[:VAR=valor,...] ..."`. `VAR` é um parâmetro do programa (a lista de cada projeto está em
  `mc/projetos.sh`; `COV_C`, `LX`, `NUCLEO`, `ANGULO`, `RHO`, `H_MALHA`, `FS_MAX`, `MALHA`, ...), e `valor` tem letras,
  números e `. _ + -`. Uma variável no ambiente vale para todas as entradas (`H_MALHA=1 mc/submeter_mc.sh V2 1`).
- **Pastas:** `$DATA/monte_carlo/<projeto>/<caso>[_trocas]/r<nnnnn>/`, com as trocas em ordem alfabética, cada uma com o
  nome da variável em minúsculas e sem `_`: `cho1:LX=40,COV_C=0.5` → `cho1_covc0.5_lx40`; `s3:ANGULO=90` →
  `s3_angulo90`. `SAIDA` troca `$DATA/monte_carlo`.
- **Determinístico:** `<caso>[_trocas]_det/r00000/` (`CAMPOS=0`, as médias), sem as trocas que só mudam os campos
  aleatórios, que ele ignora: `cho1:NUCLEO=exp` e `cho1` usam o mesmo `cho1_det`; `s3:ANGULO=90` usa `s3_det`. Cada
  pasta tem `deterministico.txt`, com o nome do seu.
- **V4, cenário s2** (k_s constante): o regime permanente é o do determinístico, com os mesmos parâmetros hidráulicos.
  Depois que ele termina, as realizações de s2 leem o `pw.csv` dele (`PW_PERMANENTE`) em vez de refazer o transiente.
- **A mesma realização tem os mesmos números aleatórios em todas as entradas** de um projeto (mesma `SEMENTE` e
  `REAL`): as diferenças entre entradas têm muito menos ruído que cada P_f sozinho (`mc/compara_mc.py`).
- **Parâmetros:** cada pasta guarda em `parametros.txt` os parâmetros com que foi rodada. Outra rodada na mesma pasta
  com outros valores (no ambiente ou nas trocas) misturaria resultados: é recusada (use outra `SAIDA`, ou uma entrada com
  as trocas). Se um projeto passar a ler uma variável que `mc/projetos.sh` não conhece, `submeter_mc.sh` avisa.
- **Retomável e seguro com vários jobs:** cada realização em andamento tem um dono (`r<nnnnn>/.ativo/dono`: máquina,
  processo, job e vez do SLURM). A de um dono vivo é pulada; a de um dono morto (job que acabou, processo que caiu,
  máquina reiniciada) volta para a fila. Um job que o SLURM dá como acabado só conta como morto depois de `INATIVO_MIN`
  minutos (padrão 30) sem mudança na pasta (um nó em `NODE_FAIL` pode estar vivo, gravando); um job recolocado na fila
  pelo SLURM (mesmo número, outra vez) retoma as próprias realizações; sem resposta do `squeue`, nenhuma trava é
  tomada. A trava é um `mkdir`, atômico também no sistema de arquivos compartilhado.
- **Falhas:** uma realização que termina com erro (sem `resultado.txt`) volta para a fila; na `MAX_FALHAS`-ésima vez
  (padrão 3) fica marcada em `falhou.txt`. Uma realização com `resultado.txt` está concluída, mesmo com `situacao`
  diferente de ok (por exemplo `sem_equilibrio_geostatico`, que é uma ruptura). **Disjuntor:** falhas seguidas demais
  numa máquina, sem nenhum sucesso (o executável ou as bibliotecas, não as realizações), fazem a instância parar de pegar
  realizações e avisar. Para tentar de novo as marcadas: `find $DATA/monte_carlo/V2_TaludesCho -name falhou.txt -delete;
  find $DATA/monte_carlo/V2_TaludesCho -name falhas -delete`.
- **Uma campanha por projeto e `SAIDA`:** `campanha.txt` guarda o último pedido, que é o que o `status_mc.sh` mostra.
  Ponha todas as entradas num pedido só (`CASOS="cho1 cho1:NUCLEO=exp"`); para acrescentar uma entrada depois, repita
  o pedido com todas elas (o que já está concluído não roda de novo), ou use outra `SAIDA` para uma campanha à parte.
- **Limite de tempo:** cada job vai até 4 dias (`TEMPO`). 15 min antes do fim, ele deixa de começar realizações e
  submete a sua continuação, que espera ele acabar e pega o que falta. Se a cota recusar a continuação, fica o aviso
  `CONTINUACAO_PENDENTE`: rode `mc/submeter_mc.sh <projeto> <nós>` quando houver vaga.
- **Avisos:** o que impediu jobs de rodar (cota, parâmetros, bibliotecas faltando num nó, disjuntor) fica em
  `$DATA/monte_carlo/<projeto>/AVISOS`, que o `status_mc.sh` mostra.
- **Compilação:** se o executável falta ou é mais velho que os fontes, `submeter_mc.sh` compila em segundo plano e
  submete ao fim (um fonte trocado durante a compilação faz compilar de novo). Recompilar com jobs `verif-*` na fila é
  recusado; para só acrescentar nós com o executável atual: `FORCAR_SEM_COMPILAR=1 mc/submeter_mc.sh V2 2`.
- **Disco:** ~25 kB por realização no V2 e ~0,1 MB no V4 (o `pw.csv` da malha de 0,5 m): a campanha padrão do V4
  (12 000 realizações) ocupa ~1,2 GB. A cota é de 512 GB em `$DATA` e em `$HOME` (`quota -s`).
- À mão, um job (de dentro de `Projects/Cluster`):
  `PROJETO=V2 CASOS="cho1" PRIMEIRA=1 ULTIMA=1000 DET=1 sbatch --job-name=verif-mc-V2 mc/job_mc.sbatch`. Fora do SLURM
  (no PC, por exemplo): `PROJETO=V2 CASOS="cho1" DET=1 SAIDA=~/runs mc/rodar_mc.sh 1 20`.

## 7. Acompanhar

```bash
./status.sh                    # compilacao, executaveis, jobs, avisos, casos completos, ultimo teste, Monte Carlo
mc/status_mc.sh V2             # o Monte Carlo do V2 em detalhe
mc/status_mc.sh --resumo V4    # e a estatistica de cada entrada (resumo_mc.sh)
squeue -u $USER
```

`mc/status_mc.sh` mostra, por entrada: o FS do determinístico, as realizações concluídas, rodando, esperando (as
interrompidas entre elas), que falharam e o total, os minutos por realização, o FS mínimo médio e P_f = P(FS_min < 1)
até agora; o fim estimado, os avisos e, no fim, o comando `rsync` para baixar os resultados.

## 8. Parar e retomar

| O quê | Parar | Retomar |
|---|---|---|
| Monte Carlo de um projeto | `touch $DATA/monte_carlo/V2_TaludesCho/PARAR` (nenhuma realização nova começa; as em andamento terminam e os jobs acabam) ou `scancel -n verif-mc-V2` | `mc/submeter_mc.sh V2 <nós>` (apaga o `PARAR`) |
| Monte Carlo de todos os projetos | `touch $DATA/monte_carlo/PARAR` | `rm $DATA/monte_carlo/PARAR` e `mc/submeter_mc.sh <projeto> <nós>` |
| Casos completos | `scancel -n verif-caso-V5` | `./submeter_caso.sh V5 <os mesmos casos>` (os concluídos são pulados) |
| Testes | `scancel -n verif-testes` | `./submeter_testes.sh` |
| Compilação em segundo plano | `kill -- -<pid>` (o pid está em `logs/compilar.pid`) | o mesmo pedido de novo |

Uma realização interrompida (pelo `scancel` ou pelo limite de tempo) volta para a fila e é refeita do início (os
programas de V2 e V4 levam um ou dois minutos por realização e não gravam checkpoint).

## 9. Baixar para o PC

Rode **no PC**, não no fermi; pode ser no meio da rodada (o `rsync` só traz o que mudou). O `mc/status_mc.sh` e o
`status.sh` mostram o comando com o caminho de `$DATA` no fermi:

```bash
rsync -av --exclude checkpoint.bin --exclude .ativo --exclude .tomada --exclude '.parar.*' --exclude /PARAR \
    --exclude /execucoes.txt --exclude /campanha.txt <usuario>@fermi.cesup.ufrgs.br:<DATA>/monte_carlo/V2_TaludesCho/ \
    ~/projects/runs/monte_carlo/V2_TaludesCho/
rsync -av --exclude .ativo --exclude /testes/ <usuario>@fermi.cesup.ufrgs.br:<DATA>/verificacao/ ~/projects/runs/verificacao/
rsync -av <usuario>@fermi.cesup.ufrgs.br:neopz-master/Projects/Cluster/logs/ ~/projects/runs/logs_cluster/
```

`$DATA` é a partição de dados de onde o manual do fermi manda rodar os jobs: temporária e sem backup. Copie os
resultados para o PC (ou para o `$HOME`) ao terminar.

## 10. Pós-processamento do Monte Carlo (no PC)

Com numpy e matplotlib, nas pastas baixadas (ou no fermi, se houver um python3 com eles):

```bash
mc/resumo_mc.sh ~/projects/runs/monte_carlo/V2_TaludesCho/cho1          # P_f, media, desvio e correlacoes
python3 mc/posprocessa_mc.py ~/projects/runs/monte_carlo/V4_TaludeLiu/s3 --mapas
python3 mc/compara_mc.py .../V2_TaludesCho/cho1 .../V2_TaludesCho/cho1_nucleoexp --rotulos "gauss,exp"
python3 mc/sensibilidade_mc.py .../V4_TaludeLiu/s3 .../V4_TaludeLiu/s3_angulo90 .../V4_TaludeLiu/s3_rho-0.5
python3 mc/plot_campos.py <r00007/campos.csv>      # os campos de uma realizacao (rodada com CAMPOS_CSV=1)
```

- `junta_mc.sh <pasta>` junta as realizações concluídas em `mc_resumo.csv` e `mc_fs.csv` (refeitos a cada chamada;
  os outros scripts o chamam sozinhos).
- `posprocessa_mc.py <pasta>` gera em `<pasta>/posproc/`: `estatisticas.csv`, a convergência de P_f (com o IC de
  Wilson) e do FS, a densidade e a distribuição acumulada do FS mínimo, as correlações de Spearman (de 1ª ordem e
  parciais) com as propriedades que variam (as médias no solo, na massa que desliza), a dispersão, o mecanismo (a massa
  que desliza e o modo de ruptura, no V4), o diagnóstico numérico e `relatorio.md`; com `--mapas`, a frequência de cada
  ponto na massa que desliza (de `mecanismo.csv`, gravado com `MEC_CSV=1`) e o p_w médio e o desvio (de `pw.csv`). O
  determinístico entra sozinho. Nos programas com um instante só do FS (V2 seco, V4 no permanente) não há as figuras por
  instante; com vários instantes (um transiente), há também a evolução no tempo.
- `compara_mc.py` compara pastas realização a realização (os mesmos campos): a diferença do FS mínimo, a correlação de
  postos, a tabela de rupturas pareada e o teste de McNemar.
- `sensibilidade_mc.py` dá as correlações de cada pasta e, com várias pastas, a tabela por caso (N, média, desvio, CoV,
  P_f, CoV(P_f)) e as densidades.
- Uma realização com `situacao` em `SITUACOES_RUPTURA` (padrão `sem_equilibrio_geostatico,fs_abaixo_de_F0`: sem
  equilíbrio nem com a resistência mínima da busca) conta como ruptura (FS = −1, fora da média); com outra situação
  diferente de ok (o cálculo não terminou, por exemplo `permanente_nao_atingido`), fica fora das estatísticas.
  `LIMITE_FS` (ou `--limite`) troca o limiar de P_f.
- Os scripts de cada projeto também leem estas pastas: `../V2_TaludesCho/cluster/pf_v2.py <pasta>` (P_f com o IC de
  Wilson contra Liao et al.) e `../V4_TaludeLiu/cluster/estatisticas_v4.py <pasta>` (β e P_f da lognormal, eqs. 12-13 de
  Liu et al.). Eles regravam `mc_resumo.csv`, que `junta_mc.sh` refaz na chamada seguinte.

## 11. Os kits dentro de cada projeto

As pastas `cluster/` de cada projeto (`compilar_v<n>.sh`, `v<n>.sbatch`, `rodar_v2.sh`, ...) foram escritas junto com
os projetos e rodam um projeto de cada vez, cada uma com a sua pasta de build. Este pacote os substitui no fermi: uma
compilação para todos os projetos, um `bin/` com as versões, a guarda contra recompilar com jobs rodando e a fila robusta
do Monte Carlo. Usar os dois compila a biblioteca duas vezes.

## 12. Como este pacote foi conferido

Sem acesso ao fermi, tudo foi testado num sandbox com falsos: `sbatch` (que roda o job em segundo plano como o SLURM,
com dependências e o aviso do fim do tempo), `squeue`, `scancel`, `scontrol`, `module`, um `cmake` falso (configura e
"compila" os alvos, com falhas simuladas) e programas falsos que seguem o contrato do Monte Carlo e o dos testes.
Foram conferidos: a compilação (do zero, a versão com e sem `+mod`, um projeto que não compila com e sem mudança na
biblioteca, as recusas com jobs, processos e `squeue` sem resposta, duas compilações ao mesmo tempo); os testes (com o
`ctest`, sem ele, com falhas, sem `[.lento]`, teste ausente, fontes mais novos); os casos completos (pastas, trocas,
listas, recusas, `REFAZER`, travas, falha do programa); a cadeia de compilação em segundo plano dos três submeter (e a
recusa de outros pedidos durante ela, o `PARAR` criado durante ela, um fonte trocado durante ela, uma compilação que
falha); e, no Monte Carlo, as pastas e os determinísticos, a ordem intercalada, `parametros.txt`, `MAX_FALHAS`, o
disjuntor, os dois `PARAR`, duas instâncias na mesma fila, a tomada da trava de donos mortos (processo, reinício, job
acabado com a pasta parada ou mudando, job recolocado na fila, `squeue` sem resposta), dois jobs na mesma fila, a
continuação antes do limite de tempo e a sua recusa, bibliotecas faltando no nó, a cota do `sbatch`, o `status_mc.sh` e
o pós-processamento. A fila também rodou com os executáveis reais de V2 (`H_MALHA=1`, `FS_MAX=1`) e de V4 (malha de 1 m,
cenário s2, com a leitura do permanente do determinístico).
