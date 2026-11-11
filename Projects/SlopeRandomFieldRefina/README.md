# SlopeRandomFieldRefina: malha refinada em volta da banda de ruptura de cada realização

Cada realização do Monte Carlo do [SlopeRandomField](../SlopeRandomField/README.md) é calculada em duas etapas de
malha: o fator na malha grossa, uma malha refinada só em volta da banda de ruptura dessa realização e o fator de novo
nela. Em vez de dividir elementos dentro do NeoPZ, o gmsh gera uma malha conforme, refinada só em volta da banda daquela
realização, sem nós pendurados. O C++ não muda: o executável é o `SlopeRandomField`, e este projeto só tem scripts.

No exemplo saturado do Ceron, no determinístico e em 8 realizações ([comparação](#comparação-com-as-malhas-fixas)), com
a malha grossa focada de 0,375 m e elementos de 0,125 m em volta da banda, o fator ficou a 0,23% em média do da malha
focada de 0,125 m (no máximo 0,49%), em 30% do tempo dela, contando a malha grossa e o gmsh. A malha focada de
0,1875 m, num tempo 4% menor, ficou a 1,58% (no máximo 2,24%).

## Método

`srf_refinado.sh` faz as três etapas de uma realização, na pasta dela:

1. **Malha grossa**, em `grossa/`: a realização na malha `MESH` (por exemplo `ceron_h5_a45_foco375.msh`, de 0,375 m
   na região das bandas), com `mecanismo.csv` só dos pontos da banda e sem `pw.csv`, `campos.csv` e VTK.
2. **Malha refinada** (`refina_banda.py`):
   - a banda são os pontos de integração com √J2(ε_p) ≥ `REFINA_LIMIAR` (padrão 0,1) do máximo, no último estado
     convergido do fator, no instante do menor fator de `fs.csv`. É a mesma definição de banda do `BANDA_LIMIAR` do
     programa;
   - os pontos são arredondados numa grade de `REFINA_H` (padrão 0,125 m): as células da banda;
   - o `.geo` novo é o da malha grossa, com a mesma geometria, os mesmos grupos físicos e os mesmos campos de tamanho,
     mais um campo: elementos de `REFINA_H` até `REFINA_FOLGA` (padrão 0,375 m) das células, crescendo até o `hMax` da
     malha grossa a `REFINA_FOLGA` + 3 m. O campo de fundo passa a ser o mínimo do antigo e deste;
   - o gmsh gera `malha_refinada.msh`: o módulo python `gmsh` do python que roda o `refina_banda.py` ou, sem ele, um
     executável `gmsh` do `PATH` que rode (`gmsh -version`).
3. **Malha refinada, de novo**: a mesma realização (mesmos `REAL`, `CAMPOS` e `SEMENTE`) na malha refinada, na pasta da
   realização. As saídas são as do programa (`resultado.txt`, `resumo.csv`, `fs.csv`, ...), e `resultado.txt` marca a
   realização como concluída, como no `rodar_mc.sh`.

O campo aleatório é o mesmo nas duas malhas:
- o `SlopeRandomField` calcula os modos de Karhunen-Loève no retângulo que contém os nós da malha e os avalia nas
  coordenadas dos pontos de integração (`KLSeparavel.h`);
- os números aleatórios de cada campo vêm de `SEMENTE`, `REAL` e do campo (c', φ', k);
- a malha refinada tem o mesmo contorno da grossa (o mesmo `.geo`), então o retângulo, os modos e os números aleatórios
  são os mesmos. Só os pontos onde o campo é avaliado mudam. No teste com o `rodar_mc.sh`, abaixo, `KL_M` e
  `KL_variancia` de `realizacao.txt` são iguais nas duas etapas;
- nas malhas refinadas do determinístico e das realizações 0, 1, 3, 5 e 7 da comparação, refeitas dos `.geo`: a mesma
  caixa dos nós da malha grossa (0 a 25 m por 0 a 10 m), os mesmos 6 cantos, os mesmos grupos físicos, os mesmos
  elementos de contorno na base e nas laterais (só a `Superficie` tem mais elementos, perto da banda) e nenhum nó fora
  dos elementos. Nas realizações, os pontos de integração que coincidem nas duas malhas (de 6 a 21, longe da banda) têm
  os mesmos c', φ' e k em `campos.csv` (`SO_CAMPOS=1`), nos 7 algarismos gravados. As médias na área (`c_med`,
  `phi_med` e `k_mult_med` de `resumo.csv`) mudam entre as etapas em até 0,9%, 0,3% e 1,8% nas realizações 0 a 7, e
  entre a `foco375` e a `foco125` em até 0,7%, 0,4% e 2,0%, porque os pontos de integração são outros; de uma
  realização para outra (as 8, na malha grossa), o desvio padrão delas é de 13,5%, 4,8% e 19,6% da média.

Sem banda na malha grossa (fator não positivo no instante do menor fator, ou o transiente não chegou até ele), não há
onde refinar: o resultado da realização é o da malha grossa, copiado para a pasta dela, com `metodo=grossa` em
`refino.csv`.

### Retomada e limpeza

- Cada etapa pula o que já está feito:
  - a 1 tem `grossa/resultado.txt`;
  - a 2 tem `refino.csv` com os mesmos `REFINA_*` e `malha_refinada.geo`;
  - a 3 continua do `checkpoint.bin` do programa, como uma realização comum.
- Um sinal no programa ou no `refina_banda.py` termina o `srf_refinado.sh` com o código do sinal (128 + sinal), sem
  mensagem do `srf_refinado.sh` (num SIGKILL ou SIGTERM só no programa, o bash grava em `run.err` a linha `Killed` ou
  `Terminated` do comando). O `rodar_mc.sh` o trata como trataria o programa: o Ctrl-C (130) e o terminal fechado
  (129) não contam como falha; um SIGKILL só no programa (137, como o do OOM killer) conta.
- Depois do sucesso, `malha_refinada.msh`, `grossa/mecanismo.csv` e um `checkpoint.bin` que tenha ficado em `grossa/`
  (de uma queda entre o `resultado.txt` e a limpeza do programa) são apagados. O `.geo` fica: se a realização for
  rodada de novo (com `resultado.txt` apagado), a malha é refeita dele. Se `refino.csv` também for apagado, a etapa 1 é
  refeita, porque a banda já não existe.
- Uma janela que fica: uma queda nos milissegundos entre o programa da etapa 3 gravar `resultado.txt` e o
  `srf_refinado.sh` regravar `refino.csv` deixa a realização concluída (o fator é o da malha refinada), mas com o
  `refino.csv` sem as colunas da etapa 3 e com `malha_refinada.msh` e `grossa/mecanismo.csv` ainda na pasta. O fator, o
  tempo e as equações dessa etapa estão em `resultado.txt` e `realizacao.txt`.
- Os `REFINA_*` de uma combinação ficam em `refino_parametros.txt`, na pasta da combinação (ao lado do
  `parametros.txt` do `rodar_mc.sh`). Outra rodada na mesma pasta com outros valores é recusada: use outra `SAIDA`.

## Arquivos

| Arquivo | O que faz |
|---|---|
| `srf_refinado.sh` | as três etapas de uma realização; é o `EXE` do `rodar_mc.sh` |
| `refina_banda.py` | a etapa 2: `refina_banda.py <rodada grossa> <malha grossa> <saída> [h] [folga] [limiar]`; `--malha <saída>` refaz o `.msh` do `.geo`; `--gmsh` mostra qual gmsh será usado (código 2 se nenhum funciona) |
| `compara.py` | a comparação de tempo e precisão abaixo: `compara.py rodar` e `compara.py tabela` |

Na pasta de cada realização, além das saídas do programa:

| Arquivo | Conteúdo |
|---|---|
| `grossa/` | a rodada na malha grossa: `resultado.txt`, `resumo.csv`, `fs.csv`, `fs_tentativas.csv`, `passos.csv`, `realizacao.txt`, `triangulos.csv`, `run.log`, `run.err` |
| `malha_refinada.geo` | o `.geo` da malha refinada (cerca de 20 kB) |
| `refino.csv` | uma linha: `real`, `campos`, `malha_grossa`, `fator_grossa`, `seg_grossa`, `neq_grossa` (equações), `celulas_banda`, `h`, `folga`, `limiar`, `triangulos` (da malha refinada), `seg_malha` (o `refina_banda.py` com o gmsh), `metodo` (`refinada` ou `grossa`), `fator`, `seg_fina`, `neq` e `seg_total` (as três etapas) |

O `seg` de `resultado.txt` e de `resumo.csv` é só o da etapa 3; o tempo das três etapas está em `seg_total`. Por isso
a coluna `min/real` do `status_mc.sh` e o tempo por realização do `resumo_mc.sh`, que leem esse `seg`, não contam a
malha grossa nem o gmsh.

`junta_mc.sh` junta os `resumo.csv` e `fs.csv` da etapa 3, como numa campanha comum. Os `refino.csv` de uma combinação
numa tabela só:

```bash
awk 'FNR > 1 || NR == 1' <SAIDA>/f50_I0_iii/r*/refino.csv > <SAIDA>/f50_I0_iii/mc_refino.csv
```

## Rodar

Precisa do executável do `SlopeRandomField` (`../SlopeRandomField/bin/SlopeRandomField`, compilado com
`../SlopeRandomField/compilar.sh`), de python3 e do gmsh: o módulo python (`python3 -m pip install --user gmsh`) ou um
executável `gmsh` no `PATH`. O módulo vem primeiro; um executável só é usado se o módulo faltar ou não carregar, e se
`gmsh -version` rodar. Para conferir qual será usado:

```bash
python3 ../SlopeRandomFieldRefina/refina_banda.py --gmsh
```

`srf_refinado.sh` entra no lugar do executável no `rodar_mc.sh`, com `EXE`; `MALHA` é a malha grossa, que precisa do
`.geo` ao lado (as malhas do Ceron em `../LearningMohrCoulomb2/cluster/malhas` têm). Tudo o mais é como no
`SlopeRandomField`: as variáveis do modelo, `COMBS`, `DET`, `PARALELO`, `SAIDA`, a retomada, `junta_mc.sh`,
`resumo_mc.sh` e `status_mc.sh`. Use uma `SAIDA` própria: o `parametros.txt` do `rodar_mc.sh` não distingue uma rodada
com o executável de uma com `srf_refinado.sh`.

Variáveis de `srf_refinado.sh`:

| Variável | Padrão | O que é |
|---|---|---|
| `SRF_BIN` | `../SlopeRandomField/bin/SlopeRandomField` | o executável |
| `REFINA_H` | 0.125 | tamanho dos elementos em volta da banda (m) |
| `REFINA_FOLGA` | 0.375 | distância da banda até onde vale `REFINA_H` (m); depois, o tamanho cresce até o `hMax` da malha grossa em 3 m |
| `REFINA_LIMIAR` | 0.1 | a banda: √J2(ε_p) ≥ limiar × máximo |
| `PYTHON` | python3 | o python de `refina_banda.py` (o do módulo `gmsh`) |

### Localmente

O exemplo saturado do Ceron, como nas fases do fermi do `SlopeRandomField`, com a malha grossa
`ceron_h5_a45_foco375` e o skyline (máquina sem MKL), o determinístico e as realizações 0 a 9:

```bash
cd Projects/SlopeRandomField
python3 -m pip install --user gmsh   # uma vez, se ainda nao houver gmsh
export EXE=$PWD/../SlopeRandomFieldRefina/srf_refinado.sh
export VARIANTE=iii MALHA=ceron_h5_a45_foco375.msh PHI_DEG=30 HL=10 HR=10 XL=0 XR=25 XTOE=15 YTOE=5 XCRISTA=10 YCRISTA=10 NIVEL_FINAL=5 T_REBAIX_H=0 FATOR=gim_arco FS_H=0.01 DT_H=0.01 SATURADO=1 SUPERFICIE_P0=1 KW=1e30 SOLVER=skyline
SAIDA=$PWD/resultados_mc_refina COMBS="f50:0" DET=1 PARALELO=2 ./rodar_mc.sh 0 9
./resumo_mc.sh resultados_mc_refina/f50_I0_iii
```

Interrompida (Ctrl-C, queda da máquina), a campanha continua com o mesmo `./rodar_mc.sh 0 9`.

### No fermi

```bash
cd ~/neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomField
```

O gmsh, uma vez, no nó de login (o `$HOME` é o mesmo nos nós de computação), e a conferência no nó de login e num nó de
computação:

```bash
python3 -m pip install --user gmsh
python3 ../SlopeRandomFieldRefina/refina_banda.py --gmsh
srun -n 1 python3 ../SlopeRandomFieldRefina/refina_banda.py --gmsh
```

- O `--gmsh` mostra o módulo e o python que o carrega, ou o executável. O `pip install --user` instala o módulo para a
  versão do python3 que o rodou (em `~/.local/lib/python3.X`): um nó de computação com outra versão do python3 não o
  encontra.
- A roda do pip do gmsh 4.15.2 é `manylinux_2_24`: precisa da glibc 2.24 ou mais nova e de um pip recente. Num
  sistema mais antigo, o pip não acha versão que instale.
- O `libgmsh.so` do pacote do pip (4.15.2) depende de bibliotecas do X e do OpenGL (libGLU, libGL, libX11, libXrender,
  libXcursor, libXfixes, libXext, libXft, libXinerama e libfontconfig). Sem uma delas o módulo não carrega, e o
  `--gmsh` mostra a mensagem do carregador, com o nome da biblioteca.
- O pip instala também `~/.local/bin/gmsh`, um script com `#!/usr/bin/env python`. Num nó sem um `python` no `PATH`
  (só `python3`), ele não roda (`/usr/bin/env: 'python': No such file or directory`). Por isso o `refina_banda.py` usa
  o módulo primeiro, e um executável só se `gmsh -version` rodar. Os dois geram o mesmo `.msh` (conferido byte a byte
  com o gmsh 4.15.2, nas malhas refinadas do determinístico e da realização 0 da comparação abaixo).
- Sem gmsh que funcione num nó, cada realização sai logo no início com código 2, antes da malha grossa, com a mensagem
  no `run.err`. Depois de `DISJUNTOR` falhas seguidas (padrão: o maior entre 2 × `PARALELO` e `MAX_FALHAS` + 2), o
  disjuntor do `rodar_mc.sh` para a instância. As realizações ficam com uma tentativa em `falhas` e voltam à fila
  quando a campanha é submetida de novo.

Com o executável já compilado (`bin/versao.txt`), as variáveis da fase e o `EXE`, exportados: o `submeter_mc.sh` os
passa aos jobs, como as variáveis do modelo.

```bash
export EXE=$PWD/../SlopeRandomFieldRefina/srf_refinado.sh
export VARIANTE=iii PRIMEIRA=0 ULTIMA=999 DET=1 MALHA=ceron_h5_a45_foco375.msh PHI_DEG=30 HL=10 HR=10 XL=0 XR=25 XTOE=15 YTOE=5 XCRISTA=10 YCRISTA=10 NIVEL_FINAL=5 T_REBAIX_H=0 FATOR=gim_arco FS_H=0.01 DT_H=0.01 SATURADO=1 SUPERFICIE_P0=1 KW=1e30
```

```bash
SAIDA=$DATA/ceron_sat_refina COMBS="f50:0" FORCAR_SEM_COMPILAR=1 ./submeter_mc.sh 4
```

```bash
SAIDA=$DATA/ceron_sat_refina ./status_mc.sh
```

- O `job_mc.sbatch` confere `bin/SlopeRandomField` e as suas bibliotecas, que é o executável que `srf_refinado.sh`
  usa.
- A `versao` de `resumo.csv` é a do executável (`bin/versao.txt`): o `rodar_mc.sh` procura `versao.txt` ao lado do
  `EXE`, e `srf_refinado.sh` a lê ao lado de `SRF_BIN`.
- O `min/real` do `status_mc.sh` é só o da etapa 3 (ver "Arquivos"); o tempo inteiro de cada realização está em
  `seg_total`, no `refino.csv` dela.
- As outras fases (trocas em `COMBS`, `NIVEL_FINAL`, ...) rodam do mesmo jeito, cada uma numa `SAIDA` própria.

## Comparação com as malhas fixas

O exemplo saturado do Ceron, com as variáveis de "Localmente" (as das fases do fermi do `SlopeRandomField`): Γ com o
GIM com arc-length em 0,01 h, no determinístico (`det`, as médias) e nas realizações 0 a 7, cada caso com cinco métodos:

- `foco375`, `foco25`, `foco1875` e `foco125`: o programa direto nas malhas `ceron_h5_a45_foco*`, com elementos de
  0,375, 0,25, 0,1875 e 0,125 m numa região fixa, a mesma nas quatro (as bandas do determinístico e das realizações 0 a
  31 do exemplo do Ceron, mais 0,75 m; ver "Malhas focadas" no README do `SlopeRandomField`);
- `refinada`: `srf_refinado.sh` com a malha grossa `foco375` e os padrões (`REFINA_H=0.125`, `REFINA_FOLGA=0.375`,
  `REFINA_LIMIAR=0.1`): em volta da banda, os elementos têm o tamanho dos da `foco125`.

A referência é a `foco125`; a diferença de um método num caso é (Γ − Γ_foco125) / Γ_foco125. Os comandos, na pasta
deste projeto:

```bash
SOLVER=skyline python3 compara.py rodar <pasta> --casos det,0-7 --paralelo 2
python3 compara.py tabela <pasta> --csv <pasta>/tabela.csv
```

### Como os tempos foram medidos

- **Tempo de parede:** o relógio do python (`time.time()`), de logo antes de o `compara.py` criar o processo do caso até
  o `wait4` que o recolhe. É o processo inteiro: carregar o programa, ler a malha, os modos de Karhunen-Loève, o estado
  inicial, o transiente, a busca do fator e a gravação das saídas. No método de duas etapas o processo é o
  `srf_refinado.sh`, e o tempo inclui a etapa 1 (a malha grossa), a etapa 2 (o `refina_banda.py` com o gmsh), a etapa 3
  (a malha refinada) e o bash.
- **Tempo de CPU:** usuário + sistema do processo e dos descendentes que ele esperou (`ru_utime + ru_stime` do
  `wait4`). Cada caso é um processo de 1 thread (`MKL_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `SOLVER=skyline`).
- **Etapas** (de `refino.csv`): a malha grossa e a refinada pelo `seg` do programa, o relógio `steady_clock` do início
  do `main` até o resumo final (sem carregar o executável e as bibliotecas); o gmsh pelo relógio do `srf_refinado.sh`
  (`date +%s.%N`) antes e depois da chamada do `refina_banda.py`: iniciar o python, ler a banda, escrever o `.geo`,
  gerar a malha e contar os triângulos. O tempo de parede é a soma das etapas mais o que fica fora desses relógios
  (carregar o programa duas vezes, o bash).
- **Máquina:** 4 núcleos (Intel Xeon a 2,10 GHz, uma thread por núcleo), 15 GB, sem MKL. O mesmo executável e a mesma
  `libpz.so` em todos os casos (`1722fa5+mod 2026-11-10_19:02:44`: o executável de `../SlopeRandomField/bin` e a
  `libpz.so` do build dele naquele momento), copiados para uma pasta à parte antes da rodada e usados com `SRF_BIN` e
  `LD_LIBRARY_PATH`, porque a biblioteca do build seria recompilada durante a rodada.
- **Fila:** os 45 casos numa fila só, 2 processos ao mesmo tempo, intercalada: para cada caso (det, 0, 1, ..., 7), os
  cinco métodos na ordem `foco125`, `refinada`, `foco1875`, `foco25` e `foco375`. Uma variação da carga da máquina ao
  longo da rodada pesa assim sobre todos os métodos. "processos juntos", na tabela, é a média no tempo de quantos casos
  da fila rodavam junto com cada um, contando ele mesmo.
- **Carga:** a máquina era compartilhada com outros trabalhos. A carga média de 1 min (`/proc/loadavg`) no início e no
  fim de cada caso está no `tempo.txt` dele (com só os 2 processos da fila, ela fica perto de 2): de 1,1 a 4,2 nos det
  e de 1,6 a 7,8 nas realizações. Amostrada a cada 30 s de 02:08 até o fim da fila (03:09), ela foi de 1,7 a 7,7, e
  acima de 2,5 em 107 das 124 amostras. Os outros processos eram compilações (de 01:54 a 02:09, de 02:15 a 02:31 e por
  1 min às 02:43), um python3 com cerca de 1,3 núcleo (de 02:27 a 02:30) e dois outros programas, cada um com até 2
  núcleos (de 02:36 a 03:05). Com mais processos prontos que núcleos, o do caso esperou: o tempo de parede passou o de
  CPU em 50% na `foco25` da realização 3 (105,5 s contra 70,3 s; carga de 7,3 a 7,7) e em no máximo 17% nos outros
  casos (a `foco1875` da realização 3). Por isso as tabelas trazem os dois tempos.
- **Duas vezes:** a rodada foi interrompida uma vez, quando o contêiner reiniciou. Os det dos cinco métodos são da
  primeira vez (01:41 a 01:50); as realizações 0 a 7, da segunda (01:56 a 03:10). Os dois casos que rodavam na
  interrupção (`foco125` e `refinada` da realização 0) recomeçaram do zero: o `compara.py` apaga a pasta de um caso sem
  `tempo.txt`, para o tempo ser o de uma execução inteira.
- **gmsh:** a `refinada` do det e da realização 0 usou o executável `gmsh` do pip (a versão anterior do
  `refina_banda.py` preferia o executável); as outras, o módulo. A malha é a mesma: o `.msh` do det e o da realização 0,
  refeitos pelos dois caminhos, são iguais byte a byte. O gmsh levou de 0,16 a 0,20 s, exceto na realização 0 (4,3 s).

### Resultados

Γ de cada caso e a diferença para a `foco125`:

| caso | foco375 | foco25 | foco1875 | refinada | foco125 |
|---|---|---|---|---|---|
| det | 1,4633 (+5,20%) | 1,4210 (+2,15%) | 1,4138 (+1,64%) | 1,3976 (+0,47%) | 1,3910 |
| 0 | 1,1224 (+5,22%) | 1,0883 (+2,02%) | 1,0845 (+1,66%) | 1,0644 (−0,21%) | 1,0667 |
| 1 | 1,4171 (+3,74%) | 1,3889 (+1,67%) | 1,3851 (+1,40%) | 1,3727 (+0,49%) | 1,3660 |
| 2 | 1,1731 (+6,15%) | 1,1360 (+2,80%) | 1,1298 (+2,24%) | 1,1029 (−0,20%) | 1,1051 |
| 3 | 1,9124 (+3,96%) | 1,8688 (+1,59%) | 1,8554 (+0,86%) | 1,8413 (+0,09%) | 1,8395 |
| 4 | 1,4702 (+6,79%) | 1,4095 (+2,38%) | 1,4036 (+1,95%) | 1,3787 (+0,14%) | 1,3767 |
| 5 | 0,9544 (+5,70%) | 0,9216 (+2,06%) | 0,9184 (+1,70%) | 0,9037 (+0,08%) | 0,9030 |
| 6 | 1,5654 (+7,23%) | 1,5067 (+3,21%) | 1,4897 (+2,04%) | 1,4642 (+0,29%) | 1,4599 |
| 7 | 1,1116 (+2,54%) | 1,0896 (+0,51%) | 1,0918 (+0,72%) | 1,0852 (+0,11%) | 1,0840 |

Tempo de parede de cada caso (s); na `refinada`, entre parênteses, as etapas: malha grossa + gmsh + malha refinada:

| caso | foco375 | foco25 | foco1875 | refinada | foco125 |
|---|---|---|---|---|---|
| det | 15,6 | 52,9 | 161,1 | 246,2 (16,0 + 0,2 + 229,8) | 551,1 |
| 0 | 12,8 | 47,8 | 175,9 | 163,5 (12,7 + 4,3 + 146,2) | 519,8 |
| 1 | 14,3 | 47,3 | 180,1 | 215,2 (14,7 + 0,2 + 200,1) | 652,4 |
| 2 | 14,3 | 56,9 | 189,4 | 167,0 (14,5 + 0,2 + 152,1) | 495,0 |
| 3 | 20,7 | 105,5 | 210,0 | 225,8 (20,4 + 0,2 + 204,9) | 981,3 |
| 4 | 16,0 | 49,5 | 171,3 | 227,6 (16,1 + 0,2 + 210,9) | 587,3 |
| 5 | 15,7 | 49,1 | 183,3 | 145,8 (13,6 + 0,2 + 131,7) | 467,9 |
| 6 | 14,5 | 51,6 | 217,9 | 235,5 (14,8 + 0,2 + 220,2) | 853,3 |
| 7 | 15,4 | 56,8 | 171,4 | 107,7 (15,2 + 0,2 + 92,1) | 601,3 |

Por método, nos 9 casos (as diferenças em valor absoluto, menos a "com sinal"; os tempos, médias por caso):

| método | triângulos | equações | dif. média | dif. máxima | dif. média com sinal | Γ < 1 | tempo de parede (s) | tempo de CPU (s) | parede / foco125 | processos juntos |
|---|---|---|---|---|---|---|---|---|---|---|
| foco375 | 720 | 3808 | 5,17% | 7,23% | +5,17% | 1 | 15,5 | 15,1 | 0,02 | 2,00 |
| foco25 | 1520 | 7898 | 2,04% | 3,21% | +2,04% | 1 | 57,5 | 53,3 | 0,09 | 1,93 |
| foco1875 | 2906 | 14960 | 1,58% | 2,24% | +1,58% | 1 | 184,5 | 179,7 | 0,29 | 2,00 |
| refinada | 2385 | 12199 | 0,23% | 0,49% | +0,14% | 1 | 192,7 | 186,0 | 0,30 | 2,00 |
| foco125 | 5426 | 27650 | 0,00% | 0,00% | +0,00% | 1 | 634,4 | 616,7 | 1,00 | 1,98 |

Na `refinada`, triângulos e equações são os da malha refinada (médias). As etapas da `refinada`, em média por caso
(a última coluna: os pontos da banda final a mais de `REFINA_FOLGA` das células da banda da malha grossa, e a maior
distância deles a essas células):

| fator na malha grossa: dif. média | células de banda | triângulos refinados | equações refinadas | malha grossa (s) | gmsh (s) | malha refinada (s) | sem banda | banda final fora da região refinada |
|---|---|---|---|---|---|---|---|---|
| 5,17% | 225 | 2385 | 12199 | 15,3 | 0,64 | 176,4 | 0 | 0 de 9 casos (no pior, 0,0% dos pontos; distância máxima 0,16 m) |

- **Precisão:** a `refinada` ficou a 0,23% da `foco125` em média e a no máximo 0,49% (realização 1), dos dois lados
  (de −0,21% a +0,49%). As malhas fixas ficaram sempre acima da `foco125`: a `foco1875` a 1,58% em média (até
  2,24%), a `foco25` a 2,04% (até 3,21%) e a `foco375` a 5,17% (até 7,23%). O Γ da etapa 1 é o da `foco375` (os 9
  valores são iguais): a etapa 3 leva a diferença de 5,17% a 0,23%. A realização 5 tem Γ < 1 nos cinco métodos; as
  outras, Γ > 1 nos cinco.
- **Tempo:** com a malha grossa e o gmsh, a `refinada` levou 192,7 s de parede por caso (186,0 s de CPU): 30% do tempo
  da `foco125` e 4% mais que a `foco1875`, tanto no de parede quanto no de CPU. Caso a caso, a razão do tempo de parede
  para o da `foco125` foi de 0,18 a 0,45, e para o da `foco1875`, de 0,63 a 1,53; a do tempo de CPU, de 0,18 a 0,42 e
  de 0,63 a 1,38. O 1,53 é o do det, cuja `refinada` rodou com a carga subindo de 1,1 a 4,2 e passou o tempo de CPU em
  11% (na repetição abaixo, 1,38). A etapa 3 é 91,6% do tempo de parede; a malha grossa, 8,0%; o gmsh, 0,3%. O tempo
  de parede passou a soma das três etapas em 0,22 a 0,29 s.
- **Malha e iterações:** a malha refinada teve de 1799 a 2624 triângulos (de 9242 a 13412 equações), 18% menos
  equações que a `foco1875` e 56% menos que a `foco125`, em média. A busca do fator nela precisou de mais iterações de
  Newton (101,7 em média, de 76 a 143) que nas malhas fixas (de 83,7 a 93,4 em média); por iteração, 1,74 s, contra
  1,97 s da `foco1875` e 6,82 s da `foco125`. O pico de memória foi de 127 a 190 MB, contra 147 a 165 MB da `foco1875`
  (com o skyline, a memória é a do perfil da matriz).
- **Bandas:** a banda de ruptura da `foco125` (√J2 ≥ 0,1 do máximo) ficou, nos 9 casos, a no máximo 0,24 m dos pontos
  do foco da `foco125` (os elementos de 0,125 m vão até 0,75 m deles). A banda final da `refinada` ficou a no máximo
  0,16 m das células da banda da malha grossa, dentro da região de elementos de 0,125 m (até 0,375 m delas). Nenhum caso
  ficou sem banda na malha grossa.
- **Repetição:** o det e a realização 7 foram medidos de novo com os cinco métodos (`compara.py rodar <pasta> --casos
  det,7`, o mesmo executável e a mesma fila de 2 processos), de 03:24 a 03:49, com a carga de 1,45 a 3,19 no
  início e no fim de cada caso. A `foco125` dos dois casos e a `refinada` da realização 7 foram rodadas outra vez
  porque um teste rodou junto com elas por 2 min; os tempos são os da segunda vez. Os fatores foram os mesmos, em
  todos os algarismos. Por caso, o tempo de CPU mudou de −1,9% a +3,6%, e o de parede de −6,7% (a `refinada` do det,
  229,7 s contra 246,2 s) a +3,9%. Nos dois casos juntos, a `refinada` levou 0,291 do tempo de parede da `foco125`
  (0,307 na primeira vez) e 0,995 do da `foco1875` (1,064); no de CPU, 0,293 (0,294) e 0,992 (0,995).

As tabelas são as de `compara.py tabela` (com `--csv`, os dados de cada caso num CSV).

## Teste com o `rodar_mc.sh`

Uma campanha pequena, com a mesma cópia do executável da comparação (`SRF_BIN` e `LD_LIBRARY_PATH`) e os scripts antes
das duas últimas mudanças (a tolerância do instante da banda em `refina_banda.py` e a limpeza de
`grossa/checkpoint.bin`, que não mudam este caso): o exemplo saturado do Ceron (as variáveis de "Localmente"), com
`COMBS="f50:0" DET=1 PARALELO=2` e `./rodar_mc.sh 0 1` (o determinístico e as realizações 0 e 1). O `rodar_mc.sh`
rodou numa sessão própria, e cada sinal foi para o grupo de processos inteiro, como num terminal. A campanha foi
interrompida duas vezes e retomada com o mesmo comando:

1. `kill -9` no grupo (como uma queda da máquina), 30 s depois de o determinístico e a realização 0 entrarem na etapa 3.
   Ficaram as travas `.ativo`, `malha_refinada.msh` e o `refino.csv` sem as colunas da etapa 3. Na retomada, os dois
   tomaram as travas dos donos mortos, pularam as etapas 1 e 2 e refizeram a etapa 3 desde o início: com um passo só
   até o instante do fator, o programa não grava checkpoint antes do fator.
2. SIGINT no grupo (como um Ctrl-C), 6 s depois de a realização 1 entrar na etapa 1, com o determinístico na etapa 3 e a
   realização 0 já concluída. Todos os processos terminaram em até 3 s, sem contar falha (nenhum `falhas`). Na
   retomada, o determinístico refez a etapa 3, e a realização 1 refez a etapa 1 (não havia `grossa/resultado.txt`) e
   seguiu até o fim.

No fim, as três tinham `resultado.txt` e `refino.csv` completo; não havia `falhas`, `falhou.txt`, `.ativo`,
`checkpoint.bin`, `malha_refinada.msh`, `grossa/mecanismo.csv` nem `AVISOS`, e os `run.err` estavam vazios;
`mc_resumo.csv` e `mc_fs.csv` do `junta_mc.sh` tinham as três; a `versao` dos `resumo.csv` era a do executável; e
`KL_M` = 2000 e `KL_variancia` = 0,99205 nas duas etapas das realizações. Os fatores (1,397625449 no determinístico,
1,064449027 e 1,372722834 nas realizações 0 e 1), os da malha grossa e os `malha_refinada.geo` são iguais aos da
comparação, rodada sem interrupção.

Uma segunda campanha, com os scripts desta versão e a mesma cópia do executável, `./rodar_mc.sh 5 7` (o determinístico
e as realizações 5 a 7, `PARALELO=2`), foi interrompida em cada etapa e retomada com o mesmo comando:

1. `kill -9` no grupo com o determinístico e a realização 5 na etapa 1: na retomada, os dois refizeram a etapa 1;
2. `kill -9` no grupo com os dois na etapa 2, com `malha_refinada.geo` e `malha_refinada.msh` gravados e sem
   `refino.csv` (um `PYTHON` que espera 25 s depois do `refina_banda.py`): na retomada, os dois pularam a etapa 1 e
   refizeram a etapa 2;
3. 20 s depois de os dois entrarem na etapa 3, SIGKILL só no programa do determinístico (como o OOM killer) e SIGINT só
   no programa da realização 5. O `rodar_mc.sh` contou uma falha no determinístico (código 137, `falhas` = 1 e a linha
   `Killed` do bash no `run.err`) e uma interrupção na realização 5 (sem `falhas`), e passou às realizações 6 e 7; com
   as duas na etapa 1, SIGINT no grupo;
4. a retomada até o fim: o determinístico e a realização 5 refizeram só a etapa 3, e a 6 e a 7, as três.

No fim, as 62 verificações deram certo: `resultado.txt` e `refino.csv` completo nas quatro, com `seg_total` igual à
soma das etapas; os fatores, os da malha grossa e os `malha_refinada.geo` iguais aos da comparação; `KL_M` e
`KL_variancia` iguais nas duas etapas das realizações; a `versao` do executável em `resumo.csv` e nas duas
`realizacao.txt`; `falhas` só no determinístico; nenhum `.ativo`, `falhou.txt`, `checkpoint.bin`,
`malha_refinada.msh`, `grossa/mecanismo.csv`, arquivo `.tmp` nem `AVISOS`; os `mc_resumo.csv` com as três realizações
e com o determinístico.

Com dois instantes de fator (`FS_H=0.01,0.02`, só a realização 7, `PARALELO=1`), o programa grava `checkpoint.bin`
depois do primeiro. Um `kill -9` no grupo durante o segundo fator da etapa 1 e outro durante o segundo fator da etapa
3: na retomada, cada etapa continuou do seu checkpoint (`retomando do checkpoint: t = 0.01 h` em `grossa/run.log` e em
`run.log`; `retomadas` = 1 e 2 passos em `grossa/resumo.csv` e em `resumo.csv`), com a `malha_refinada.geo` e o fator
(1,085179691) da comparação.

Com o executável atual de `../SlopeRandomField/bin`, sem `SRF_BIN` nem `LD_LIBRARY_PATH` (a versão
`ddfe3bb+mod 2026-11-11_02:09:10`, com a `libpz.so` do build com as mudanças do `TPZKarhunenLoeveMat` e do
`TPZInterpolationSpace`), `./rodar_mc.sh 7 7` (`DET=0`, `PARALELO=1`) deu os fatores da comparação (1,111594024 na
malha grossa e 1,085179691 na refinada), a mesma `malha_refinada.geo`, o mesmo `fs.csv` (fora a coluna `seg`) e a
`versao` desse executável em `resumo.csv`.

Sem gmsh (um python sem o módulo e nenhum executável `gmsh` no `PATH`), com o executável falso descrito abaixo,
`PARALELO=2` e as realizações 0 a 9: cada realização saiu com código 2 antes da malha grossa, com a mensagem no
`run.err`, e o disjuntor parou a instância depois de 5 falhas seguidas (uma sexta já tinha começado). As seis ficaram
com uma tentativa em `falhas`. Sem gmsh, com o executável congelado, `MAX_FALHAS=2`, `DISJUNTOR=100` e as realizações 0
a 2: cada uma saiu com código 2 na primeira passada, voltou à fila na segunda, saiu de novo com código 2 e ficou
marcada em `falhou.txt` (`codigo 2 em 2 tentativas`).

Os outros ramos de `srf_refinado.sh` foram conferidos com um executável falso, que na malha grossa copia as saídas de
uma rodada real (a `foco375` do determinístico da comparação) e na refinada escreve um resultado fixo; as 36
verificações deram certo:

- a realização com banda: a malha refinada com as mesmas 271 células e 2580 triângulos da comparação, `refino.csv`,
  `refino_parametros.txt`, `seg_total` igual à soma das etapas e a limpeza do fim;
- rodada de novo depois de apagar `resultado.txt`: só a etapa 3, com a malha refeita do `.geo`; depois de apagar também
  `refino.csv`: as três etapas, porque a banda tinha sido apagada;
- sem banda na malha grossa: o resultado da malha grossa, com `metodo=grossa`, sem a etapa 3;
- outro `REFINA_H` na mesma combinação: recusado (código 2), sem rodar nada;
- um sinal na etapa 1: código 130, sem mensagem de erro;
- a malha grossa sem `.geo` e `REFINA_LIMIAR` fora de (0, 1]: código 2;
- `SO_CAMPOS=1`: o executável direto, na malha grossa, sem `grossa/` nem `refino_parametros.txt`;
- sem o módulo `gmsh`, com um executável `gmsh` que roda: a mesma `malha_refinada.geo` e os mesmos 2580 triângulos;
- sem o módulo `gmsh`, com o executável do pip num `PATH` sem `python`: código 2 antes da malha grossa, com
  `nao roda (codigo 127)` no `run.err`.

## Estudo dos níveis de refinamento e relatório

`niveis.py` roda o estudo dos níveis (o tamanho `h` na banda e a malha inicial da etapa 1) no exemplo saturado do
Ceron, nos mesmos 9 casos da comparação (det e realizações 0 a 7), com um executável congelado e o tempo medido como em
`compara.py`. A rodada do relatório usou as malhas iniciais `foco375` e `ceron_h5_a45_g15` (126 triângulos, nova em
`../LearningMohrCoulomb2/cluster/malhas`) e os níveis 0,25, 0,1875 e 0,125 m, com 3 processos ao mesmo tempo:

```bash
python3 niveis.py rodar <pasta> --paralelo 3 --sem-extra --configs foco375/um_passo_h0.25,foco375/um_passo_h0.1875,g15/um_passo_h0.25,g15/um_passo_h0.1875,g15/um_passo_h0.125
```

O relatório (`relatorio/relatorio.pdf`, fontes em `relatorio/`) junta a comparação, os níveis e o Karhunen-Loève
(formulações, meio cilindro de Mika et al. 2021 e talude). Recomendação: a `foco375` como malha inicial e um passo até
h = 0,125 m (a configuração padrão de `srf_refinado.sh`). Para refazer as tabelas e figuras: `python3
relatorio/analise_niveis.py`, `python3 relatorio/figuras_malhas.py todas`, `python3 relatorio/tabelas_comparacao.py` e
os scripts de `../KLCilindroMika/figuras/` (leem os dados das rodadas, que ficaram fora do repositório; os valores das
figuras e tabelas estão em `relatorio/dados/`). Compilar com `pdflatex relatorio.tex` (duas vezes) dentro de `relatorio/`.
