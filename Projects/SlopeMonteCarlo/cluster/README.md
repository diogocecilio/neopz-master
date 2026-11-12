# Seções 6.1 a 6.5 de Vargas Ceron et al. (2025) no fermi

O `cluster.cpp` (alvo `SlopeMonteCarloCluster`) roda o Monte Carlo do rebaixamento do Ceron com o mesmo modelo do
`main.cpp` deste projeto:
- malha `ceron_h5_a45_g15`, com refinamento adaptativo em h (nunif = 0, nref = 4, pmax = 2);
- percolação permanente com o k_v da realização;
- Γ pelo GIM com arc-length, no material não saturado;
- c', φ' e k_v lognormais, gerados por Karhunen-Loève com M = 2000.

A realização n tem os mesmos números aleatórios em todos os casos.

O padrão é a campanha de 10 000 realizações por caso com 4 níveis de refinamento adaptativo. A primeira campanha, de
1000 realizações com 3 níveis, fica na pasta dela (abaixo, "A campanha de 3 níveis").

## Casos

Os 26 casos das Tabelas 3 a 6 do artigo (`bin/SlopeMonteCarloCluster casos` lista):

| Seção | Casos | Pasta |
|---|---|---|
| 6.1 | referência (Tabela 2) e Spearman | `<prefixo>_fase1/f50_I0_iii` |
| 6.2 | CoV(k_v) 0; 0,75; 0,9; 1 · CoV(c') 0,1; 0,5; 0,7 · CoV(φ') 0,05; 0,15; 0,2 | `<prefixo>_fase2/f50_I0_iii_covk0`, ... |
| 6.3 | s = 1,5; 2; 5; 10; 20; 400 (L_x = 20s m, L_y = 2s m) | `<prefixo>_fase2/f50_I0_iii_lx30_ly3`, ... |
| 6.4 | α = k_h/k_v = 2, 3, 4, 5 | `<prefixo>_fase2/f50_I0_iii_alfa2`, ... |
| 6.5 | h_w/H = 0,5 a 0,9 (nível final 7,5 a 5,5 m) | `<prefixo>_fase3/f50_I0_iii_nivelfinal7.5`, ... |

- Os casos com determinístico próprio (referência, α e h_w) têm também a pasta `<caso>_det/r00000`.
- As pastas e os arquivos são os da campanha `ceron_sat` do SlopeRandomField, para que
  `../../SlopeRandomField/comparacao_ceron/dados_sat.py` e `relatorio_sat.py` leiam os resultados sem mudanças. O nome
  `f50_I0_iii` vem daquele programa (solo f50, sem chuva, variante iii) e aqui só identifica as pastas. O prefixo é
  `ceron_mc` (`SMC_PREFIXO` muda).
- No rebaixamento parcial (6.5), a água do reservatório carrega a face e o terreno do pé abaixo do nível. A carga cresce
  com o λ do GIM, como as poropressões.

## Rodar no fermi

```bash
cd ~/neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeMonteCarlo/cluster
source ../../LearningMohrCoulomb2/cluster/modulos.sh
./compilar.sh                     # no nó de login: a biblioteca inteira leva 15 a 30 min na primeira vez
./submeter.sh 4                   # 4 nós, 12 processos cada; volta em segundos, pode desconectar
./status.sh                       # jobs na fila, últimas linhas dos logs e a tabela até agora
```

- **Fila de tarefas:** primeiro os determinísticos, depois as realizações `PRIMEIRA..ULTIMA` (padrão 0..9999) de todos
  os casos, intercaladas (a 0 de todos, depois a 1...). Cada um dos `NOS × PARALELO` processos pega uma tarefa a cada
  `NOS × PARALELO`. Parada no meio, a campanha tem quase o mesmo número de realizações em todos os casos.
- **Níveis de refinamento:** `NREF` (padrão 4) é o número de níveis adaptativos; o Γ da realização é o do último.
  Cada pasta de saída é de um `NREF` só: a que já tem resultados de outro é recusada, na submissão e no início de cada
  processo, para as realizações de malhas diferentes não entrarem na mesma estatística.
- **Saída:** os resultados vão para `$DATA/ceron_mc_h4` (`ceron_mc_h<NREF>`), na partição de dados do fermi, que não tem
  backup (copie para o PC ao terminar). Sem `$DATA`, vão para `./resultados_h4`.
- **Variáveis de `submeter.sh`:** `NREF`, `SAIDA`, `PRIMEIRA`, `ULTIMA`, `CASOS` (por exemplo `CASOS=fase1` ou
  `CASOS=ref,hw=0.9`), `PARALELO` (padrão 12, os núcleos físicos do nó), `TEMPO` (padrão 4 dias, o máximo da partição
  `short`), `PARTICAO` e `CONTINUAR=0` (sem continuação no limite de tempo).
- **Limite de tempo:** 15 min antes do limite, os processos terminam a tarefa corrente e param, e o job submete a sua
  continuação: o mesmo nó da divisão, que espera ele acabar e faz o que falta. A campanha segue assim até a fila
  acabar, sem nova submissão.
  - Se o SLURM recusar a continuação (cota de jobs), o `status.sh` avisa. Rode `./submeter.sh <nós>` quando os jobs
    acabarem.
  - O job avisado com menos de 15 min de conta (um `TEMPO` curto demais) não continua.
- **Retomar:** submeter de novo, depois que os jobs acabarem, faz só o que falta. As realizações feitas estão em
  `partes/*.csv` e são puladas.
- **Mais realizações:** `ULTIMA=19999 ./submeter.sh 4` acrescenta as realizações 10 000 a 19 999. O artigo usou blocos
  de 10 000 até CoV(P_f) < 5 %.
- **Parar:** `touch $DATA/ceron_mc_h4/PARAR` faz os processos não começarem tarefas novas. `submeter.sh` apaga o
  arquivo. `scancel -n slopeMC` cancela os jobs, os que esperam também.
- **Tarefa com erro:** a que termina numa exceção do NeoPZ (um `DebugStop`) fica sem resultado, com uma linha `ERRO` no
  log do processo. A fila segue, a próxima submissão tenta de novo, e o `status.sh` mostra quantas são.
- **Não submeta com jobs `slopeMC` na fila:** a divisão das tarefas depende de `NOS` e `PARALELO`, e duas submissões
  ao mesmo tempo repetiriam tarefas. `submeter.sh` recusa. O `compilar.sh` também: os jobs usam o executável e a
  `libpz.so` que ele trocaria.

**Custo.** As 260 010 tarefas (26 casos × 10 000, mais 10 determinísticos) com 4 níveis:
- Na campanha de 3 níveis, o fermi fez 25 365 tarefas em 6 h 35 min com 4 nós × 12 processos: 45 s por tarefa.
- O quarto nível custa 2,1 vezes o tempo dos três primeiros juntos (52 realizações, as 0 e 1 dos 26 casos, medidas com
  3 e com 4 níveis na mesma máquina: 17,8 s e 37,3 s em média; de 1,5 a 3,5 vezes por realização). A malha do último
  nível passa de 900 para 1540 triângulos em média (até 2400).
- No fermi isso dá cerca de 95 s por tarefa e 6900 h de processo: perto de 6 dias com 4 nós (um limite de tempo e uma
  continuação) e de 3 dias com 8.
- `CASOS=fase1 ./submeter.sh 4` roda só a referência (Seção 6.1): 10 001 tarefas, perto de 6 h com 4 nós. Os outros
  casos entram depois, com `./submeter.sh 4`, na mesma pasta.
- Cada processo usa cerca de 150 MB.

## A campanha de 3 níveis

A primeira campanha (1000 realizações por caso, `NREF=3`) está em `$DATA/ceron_mc`, e o executável novo a lê e a
continua:

```bash
SAIDA=$DATA/ceron_mc ./status.sh                    # a tabela (a configuração do cabeçalho vem de parametros.txt)
NREF=3 SAIDA=$DATA/ceron_mc ./submeter.sh 4         # faz o que falta; com ULTIMA maior, acrescenta realizações
```

Sem `NREF=3`, o `submeter.sh` recusa aquela pasta.

## Arquivos

Em cada caso:

| Arquivo | Conteúdo |
|---|---|
| `partes/p<k>de<K>.csv` | as realizações do processo k, uma linha cada, gravada ao fim da realização |
| `mc_resumo.csv` | todas as realizações, em ordem (o `junta` do fim de cada job; `status.sh` também o refaz) |
| `parametros.txt` | malha, refinamento, estatística dos campos, α e h_w do caso |

Colunas de `mc_resumo.csv`:

| Coluna | Conteúdo |
|---|---|
| `real` | número da realização |
| `cov_c`, `cov_phi`, `cov_k`, `Lx`, `Ly`, `alfa`, `hw_H`, `nivel_final` | parâmetros do caso |
| `KL_M`, `KL_variancia` | termos do KL e fração da variância retida |
| `c_*`, `phi_*`, `k_mult_*` | média, desvio, mínimo e máximo dos campos no solo, ponderados pela área |
| `FS_min` | o Γ da realização no último nível (nome do programa antigo, lido pelos scripts) |
| `gamma_nivel0..4` | o Γ de cada nível de refinamento (até `NREF`) |
| `banda_*` | a banda de ruptura (√J2(ε_p) ≥ 0,1 do máximo): área, centroide e c', φ', k médios |
| `massa_*` | a massa que desliza (\|Δu\| do último incremento do GIM ≥ 0,5 do máximo): as mesmas grandezas |
| `triangulos`, `equacoes`, `passos`, `fim` | malha do último nível e caminho do fator |
| `situacao`, `seg`, `versao` | situação, tempo e commit do executável |

As variáveis do Spearman da Seção 6.1 são `banda_c`, `banda_phi_deg` e `massa_k_mult`.

Em `<caso>_det/r00000/`:
- `resumo.csv` (as mesmas colunas);
- `mecanismo.csv` (os pontos com √J2(ε_p) ≥ 0,01 do máximo);
- `resultado.txt`;
- `gim_unsat_nivel<k>.vtk` (cada nível, para o ParaView).

## Relatório

Traga `$DATA/ceron_mc_h4` para o PC e rode os scripts do SlopeRandomField:

```bash
rsync -av <usuario>@fermi:<DATA>/ceron_mc_h4/ ceron_mc_h4/
cd ../../SlopeRandomField/comparacao_ceron
python3 -I dados_sat.py <pasta com ceron_mc_fase*> ceron_mc estatisticas.json
python3 -I relatorio_sat.py <pasta com ceron_mc_fase*> ceron_mc - <saida> "g15, 4 níveis h"
```

`bin/SlopeMonteCarloCluster resumo <pasta>` dá a mesma tabela, sem o Spearman.

## Teste local

Rodado antes da submissão, com 3 e com 4 níveis:
- **Determinísticos:**

  | caso | 3 níveis | 4 níveis | artigo | diferença (4 níveis) |
  |---|---|---|---|---|
  | referência | 1,4032 | 1,3850 | 1,336 | +3,7 % |
  | α = 2 | 1,5621 | 1,5400 | 1,533 | +0,5 % |
  | α = 3 | 1,6592 | 1,6361 | 1,674 | −2,3 % |
  | α = 4 | 1,7262 | 1,7034 | 1,783 | −4,5 % |
  | α = 5 | 1,7776 | 1,7531 | 1,872 | −6,4 % |
  | h_w/H = 0,5 | 1,7392 | 1,7203 | 1,671 | +2,9 % |
  | h_w/H = 0,6 | 1,5591 | 1,5427 | 1,494 | +3,3 % |
  | h_w/H = 0,7 | 1,4467 | 1,4298 | 1,383 | +3,4 % |
  | h_w/H = 0,8 | 1,3853 | 1,3697 | 1,322 | +3,6 % |
  | h_w/H = 0,9 | 1,3707 | 1,3542 | 1,307 | +3,6 % |

  O mínimo do artigo em h_w/H ≈ 0,9, abaixo do rebaixamento total, aparece.
- **Níveis da referência:** Γ = 1,7444; 1,5559; 1,4488; 1,4032; 1,3850 nos níveis 0 a 4, com 126, 237, 431, 811 e
  1314 triângulos. As diferenças entre níveis caem na razão de 0,4, e a sequência extrapola para perto de 1,373: o Γ
  do método do Ceron com a mesma p_w neste domínio (1,374, no README do SlopeRandomField). O viés da malha cai então de
  2,1 % no nível 3 para 0,8 % no nível 4. O resto da diferença para o artigo (1,336) é do contorno distante da
  percolação, e não muda com a malha.
- **Realizações:** nas 52 realizações do teste de custo, os níveis 0 a 3 são idênticos com `NREF=3` e com `NREF=4`, e o
  quarto nível baixa Γ em 1,1 % em média (de 0,6 % a 1,4 %). Todas terminaram no patamar, sem erro.
- **Contra o fermi:** os 10 determinísticos de 3 níveis e cinco realizações da campanha de 3 níveis (a 970 de s = 5, a
  974 de CoV(c') = 0,5, a 985 e a 998 de CoV(φ') = 0,15 e a 966 de h_w/H = 0,8) deram aqui, no nível 3, os Γ do fermi
  com as quatro casas (1,4559; 2,3401; 0,9137; 1,0997; 1,1561).
- **Referência com 3 níveis:** o Monte Carlo de 1000 realizações (`main.cpp`, na estação) deu P_f = 10,3 %, contra
  11,5 % do artigo, com média lognormal de 1,398 (artigo 1,353).
- **Scripts:** `submeter.sh`, `job.sbatch` e `status.sh` foram rodados com um SLURM simulado (`sbatch` e `squeue`
  falsos), incluindo o aviso do limite de tempo, a continuação e a recusa dela. `dados_sat.py` leu a saída.
