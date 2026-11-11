# SlopeRandomFieldRefina: malha refinada em volta da banda de ruptura de cada realização

Cada realização do Monte Carlo do [SlopeRandomField](../SlopeRandomField/README.md) é calculada em duas etapas de
malha: o fator na malha grossa, uma malha refinada só em volta da banda de ruptura dessa realização e o fator de novo
nela. Em vez de dividir elementos dentro do NeoPZ, o gmsh gera uma malha conforme, refinada só em volta da banda daquela
realização, sem nós pendurados. O C++ não muda: o executável é o `SlopeRandomField`, e este projeto só tem scripts.

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
   - o gmsh gera `malha_refinada.msh` (o executável `gmsh` do `PATH` ou, sem ele, o módulo python `gmsh`).
3. **Malha refinada, de novo**: a mesma realização (mesmos `REAL`, `CAMPOS` e `SEMENTE`) na malha refinada, na pasta da
   realização. As saídas são as do programa (`resultado.txt`, `resumo.csv`, `fs.csv`, ...), e `resultado.txt` marca a
   realização como concluída, como no `rodar_mc.sh`.

O campo aleatório é o mesmo nas duas malhas:
- o `SlopeRandomField` calcula os modos de Karhunen-Loève no retângulo que contém os nós da malha e os avalia nas
  coordenadas dos pontos de integração (`KLSeparavel.h`);
- os números aleatórios de cada campo vêm de `SEMENTE`, `REAL` e do campo (c', φ', k);
- a malha refinada tem o mesmo contorno da grossa (o mesmo `.geo`), então o retângulo, os modos e os números aleatórios
  são os mesmos. Só os pontos onde o campo é avaliado mudam. No teste abaixo, `KL_M` e `KL_variancia` de
  `realizacao.txt` são iguais nas duas etapas.

Sem banda na malha grossa (fator não positivo no instante do menor fator, ou o transiente não chegou até ele), não há
onde refinar: o resultado da realização é o da malha grossa, copiado para a pasta dela, com `metodo=grossa` em
`refino.csv`.

### Retomada e limpeza

- Cada etapa pula o que já está feito:
  - a 1 tem `grossa/resultado.txt`;
  - a 2 tem `refino.csv` com os mesmos `REFINA_*` e `malha_refinada.geo`;
  - a 3 continua do `checkpoint.bin` do programa, como uma realização comum.
- Um sinal (Ctrl-C, `scancel`, `kill`) sai com o código do programa (128 + sinal), que o `rodar_mc.sh` não conta como
  falha.
- Depois do sucesso, `malha_refinada.msh` e `grossa/mecanismo.csv` são apagados. O `.geo` fica: se a realização for
  rodada de novo (com `resultado.txt` apagado), a malha é refeita dele.
- Os `REFINA_*` de uma combinação ficam em `refino_parametros.txt`, na pasta da combinação (ao lado do
  `parametros.txt` do `rodar_mc.sh`). Outra rodada na mesma pasta com outros valores é recusada: use outra `SAIDA`.

## Arquivos

| Arquivo | O que faz |
|---|---|
| `srf_refinado.sh` | as três etapas de uma realização; é o `EXE` do `rodar_mc.sh` |
| `refina_banda.py` | a etapa 2: `refina_banda.py <rodada grossa> <malha grossa> <saída> [h] [folga] [limiar]`; `--malha <saída>` refaz o `.msh` do `.geo`; `--gmsh` mostra qual gmsh será usado |
| `compara.py` | a comparação de tempo e precisão abaixo: `compara.py rodar` e `compara.py tabela` |

Na pasta de cada realização, além das saídas do programa:

| Arquivo | Conteúdo |
|---|---|
| `grossa/` | a rodada na malha grossa: `resultado.txt`, `resumo.csv`, `fs.csv`, `fs_tentativas.csv`, `passos.csv`, `realizacao.txt`, `triangulos.csv`, `run.log`, `run.err` |
| `malha_refinada.geo` | o `.geo` da malha refinada (cerca de 20 kB) |
| `refino.csv` | uma linha: `real`, `campos`, `malha_grossa`, `fator_grossa`, `seg_grossa`, `neq_grossa` (equações), `celulas_banda`, `h`, `folga`, `limiar`, `triangulos` (da malha refinada), `seg_malha` (o `refina_banda.py` com o gmsh), `metodo` (`refinada` ou `grossa`), `fator`, `seg_fina`, `neq` e `seg_total` (as três etapas) |

O `seg` de `resultado.txt` e de `resumo.csv` é só o da etapa 3; o tempo das três etapas está em `seg_total`. Por isso
a coluna `min/real` do `status_mc.sh`, que lê o `seg` de `resultado.txt`, não conta a malha grossa nem o gmsh.

`junta_mc.sh` junta os `resumo.csv` e `fs.csv` da etapa 3, como numa campanha comum. Os `refino.csv` de uma combinação
numa tabela só:

```bash
awk 'FNR > 1 || NR == 1' <SAIDA>/f50_I0_iii/r*/refino.csv > <SAIDA>/f50_I0_iii/mc_refino.csv
```

## Rodar

Precisa do executável do `SlopeRandomField` (`../SlopeRandomField/bin/SlopeRandomField`, compilado com
`../SlopeRandomField/compilar.sh`), de python3 e do gmsh: o executável `gmsh` no `PATH` ou o módulo python
(`python3 -m pip install --user gmsh`). Para conferir qual será usado:

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
| `PYTHON` | python3 | o python de `refina_banda.py` |

### Localmente

O exemplo saturado do Ceron, como nas fases do fermi do `SlopeRandomField`, com a malha grossa
`ceron_h5_a45_foco375` e o skyline (máquina sem MKL), o determinístico e as realizações 0 a 9:

```bash
cd Projects/SlopeRandomField
export EXE=$PWD/../SlopeRandomFieldRefina/srf_refinado.sh
export VARIANTE=iii MALHA=ceron_h5_a45_foco375.msh PHI_DEG=30 HL=10 HR=10 XL=0 XR=25 XTOE=15 YTOE=5 XCRISTA=10 YCRISTA=10 NIVEL_FINAL=5 T_REBAIX_H=0 FATOR=gim_arco FS_H=0.01 DT_H=0.01 SATURADO=1 SUPERFICIE_P0=1 KW=1e30 SOLVER=skyline
SAIDA=$PWD/resultados_mc_refina COMBS="f50:0" DET=1 PARALELO=2 ./rodar_mc.sh 0 9
./resumo_mc.sh resultados_mc_refina/f50_I0_iii
```

### No fermi

```bash
cd ~/neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomField
```

O gmsh, uma vez, no nó de login (o `$HOME` é o mesmo nos nós de computação), e a conferência:

```bash
python3 -m pip install --user gmsh && python3 ../SlopeRandomFieldRefina/refina_banda.py --gmsh
```

- O `libgmsh.so` do pacote do pip (4.15.2) depende de bibliotecas do X e do OpenGL (libGLU, libGL, libX11, libXrender,
  libXcursor, libXfixes, libXext, libXft, libXinerama, libfontconfig e libfreetype). Se um nó de computação não as
  tiver, `refina_banda.py --gmsh` mostra a que falta.
- Para conferir num nó de computação: `srun -n 1 python3 ../SlopeRandomFieldRefina/refina_banda.py --gmsh`.
- Sem o gmsh num nó, cada realização sai logo no início com código 2, antes da malha grossa, com a mensagem no
  `run.err`. Depois de 2 × `PARALELO` falhas seguidas, o disjuntor do `rodar_mc.sh` para a instância.
- Um executável `gmsh` no `PATH` dos jobs tem preferência sobre o módulo.

Com o executável já compilado (`bin/versao.txt`), as variáveis da fase e o `EXE`, que o `submeter_mc.sh` passa aos jobs:

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
- As outras fases (trocas em `COMBS`, `NIVEL_FINAL`, ...) rodam do mesmo jeito, cada uma numa `SAIDA` própria.

## Teste com o `rodar_mc.sh`

Uma campanha pequena, numa máquina sem MKL: o exemplo saturado do Ceron (as variáveis de "Localmente"), com
`COMBS="f50:0" DET=1 PARALELO=2` e `./rodar_mc.sh 0 1` (o determinístico e as realizações 0 e 1), interrompida duas
vezes e retomada com o mesmo comando:

1. `kill -9` em todo o grupo de processos (como uma queda da máquina), com o determinístico e a realização 0 já na
   etapa 3. Na retomada, os dois pularam as etapas 1 e 2 (as travas `.ativo` dos donos mortos foram tomadas) e
   refizeram a etapa 3 desde o início: com um passo só até o instante do fator, o programa não grava checkpoint antes
   do fator.
2. SIGINT em todo o grupo (como um Ctrl-C), com a realização 1 na etapa 3. Todos os processos terminaram, sem contar
   falha (nenhum `falhas`). Na retomada, ela pulou as etapas 1 e 2 e terminou.

No fim: as três com `resultado.txt` e `refino.csv` completo; nenhum `falhas`, `falhou.txt`, `.ativo`, `checkpoint.bin`,
`malha_refinada.msh` ou `grossa/mecanismo.csv`; `mc_resumo.csv` e `mc_fs.csv` do `junta_mc.sh` com as três; `versao`
do `resumo.csv` igual à de `bin/versao.txt`; `KL_M` = 2000 e `KL_variancia` = 0,99205 nas duas etapas. Os fatores
(1,397625449 no determinístico, 1,064449027 e 1,372722834 nas realizações 0 e 1) são iguais aos das rodadas sem
interrupção da comparação abaixo.

Os outros ramos de `srf_refinado.sh` foram conferidos com um executável falso, que copia as saídas de uma rodada real
na malha grossa: sem banda (o resultado da malha grossa, com `metodo=grossa`); a realização rodada de novo depois de
apagar `resultado.txt` (a malha refeita do `.geo`) ou também `refino.csv` (a etapa 1 refeita, porque a banda tinha
sido apagada); outro `REFINA_H` na mesma combinação (recusado, código 2); um sinal na etapa 1 (código 130); a malha
grossa sem `.geo` e `REFINA_LIMIAR` fora de (0, 1] (código 2); `SO_CAMPOS=1` (o executável direto, na malha grossa).
