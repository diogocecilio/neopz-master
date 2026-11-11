# Seções 6.1 a 6.5 de Vargas Ceron et al. (2025) no fermi

O `cluster.cpp` (alvo `SlopeMonteCarloCluster`) roda o Monte Carlo do rebaixamento do Ceron com o mesmo modelo do
`main.cpp` deste projeto:
- malha `ceron_h5_a45_g15`, com refinamento adaptativo em h (nunif = 0, nref = 3, pmax = 2);
- percolação permanente com o k_v da realização;
- Γ pelo GIM com arc-length, no material não saturado;
- c', φ' e k_v lognormais, gerados por Karhunen-Loève com M = 2000.

A realização n tem os mesmos números aleatórios em todos os casos.

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

- **Fila de tarefas:** primeiro os determinísticos, depois as realizações `PRIMEIRA..ULTIMA` (padrão 0..999) de todos
  os casos, intercaladas (a 0 de todos, depois a 1...). Cada um dos `NOS × PARALELO` processos pega uma tarefa a cada
  `NOS × PARALELO`.
- **Saída:** os resultados vão para `$DATA/ceron_mc`, a partição de dados do fermi, que não tem backup (copie para o
  PC ao terminar). Sem `$DATA`, vão para `./resultados`.
- **Variáveis de `submeter.sh`:** `SAIDA`, `PRIMEIRA`, `ULTIMA`, `CASOS` (por exemplo `CASOS=fase1` ou
  `CASOS=ref,hw=0.9`), `PARALELO` (padrão 12, os núcleos físicos do nó), `TEMPO` (padrão 2 dias) e `PARTICAO`.
- **Retomar:** submeter de novo, depois que os jobs acabarem, faz só o que falta. As realizações feitas estão em
  `partes/*.csv` e são puladas.
- **Mais realizações:** `ULTIMA=1999 ./submeter.sh 4` acrescenta as realizações 1000 a 1999. O artigo usou blocos de
  10 000 até CoV(P_f) < 5 %.
- **Parar:** `touch $DATA/ceron_mc/PARAR` faz os processos não começarem tarefas novas. `submeter.sh` apaga o arquivo.
- **Não submeta com jobs `slopeMC` na fila:** a divisão das tarefas depende de `NOS` e `PARALELO`, e duas submissões
  ao mesmo tempo repetiriam tarefas. `submeter.sh` recusa; `scancel -n slopeMC` cancela os jobs.

**Custo.** Na estação (i7-10700K, 1 thread por processo), cada realização leva ~19 s. As 26 010 tarefas (26 casos ×
1000, mais 10 determinísticos) somam ~140 h de CPU. No fermi, conte com até o dobro por núcleo: 4 nós × 12 processos
devem levar de 3 a 6 h.

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
| `gamma_nivel0..3` | o Γ de cada nível de refinamento |
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

Traga `$DATA/ceron_mc` para o PC e rode os scripts do SlopeRandomField:

```bash
rsync -av <usuario>@fermi:<DATA>/ceron_mc/ ceron_mc/
cd ../../SlopeRandomField/comparacao_ceron
python3 -I dados_sat.py <pasta com ceron_mc_fase*> ceron_mc estatisticas.json
python3 -I relatorio_sat.py <pasta com ceron_mc_fase*> ceron_mc - <saida> "g15, 3 níveis h"
```

`bin/SlopeMonteCarloCluster resumo <pasta>` dá a mesma tabela, sem o Spearman.

## Teste local

Rodado antes da submissão:
- **Determinísticos:**

  | caso | este | artigo | diferença |
  |---|---|---|---|
  | referência | 1,4032 | 1,336 | +5,0 % |
  | α = 2 | 1,5621 | 1,533 | +1,9 % |
  | h_w/H = 0,5 | 1,7392 | 1,671 | +4,1 % |
  | h_w/H = 0,9 | 1,3707 | 1,307 | +4,9 % |

  O mínimo do artigo em h_w/H ≈ 0,9, abaixo do rebaixamento total, aparece. Os +5 % são o viés da malha no nível 3.
- **Referência:** o Monte Carlo de 1000 realizações deste modelo (`main.cpp`, na estação) deu P_f = 10,3 %, contra
  11,5 % do artigo, com média lognormal de 1,398 (artigo 1,353).
- **Scripts:** `job.sbatch` foi rodado com as variáveis do SLURM simuladas. `dados_sat.py` e `relatorio_sat.py` geraram
  o relatório a partir da saída.
