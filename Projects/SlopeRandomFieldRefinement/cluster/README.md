# Monte Carlo da chuva no fermi

O `montecarlo.cpp` (alvo `SlopeRainMonteCarlo`) roda o Monte Carlo do talude de Rahardjo et al. (2007) sob chuva:
- talude da Fig. 2 do Rahardjo, com H = 10 m a 45°, na `TriGMesh`;
- solo f50: k_s = 1e-5 m/s, Fredlund-Xing com a = 50 kPa, c' = 10 kPa, φ' = φ^b = 26°;
- transiente da chuva com esqueleto rígido e ar passivo, como no SEEP/W do artigo;
- FS pelo SRM, com refinamento adaptativo em h refeito do zero em cada instante (`TriGMesh(1)` + 4 níveis);
- c', φ' e k_s lognormais, gerados por Karhunen-Loève com M = 2000 e as estatísticas da Tabela 2 de Vargas Ceron et al.
  (2025).

A realização n tem os mesmos números aleatórios em todos os casos. Os casos são o análogo das Seções 6.1 a 6.5 do Ceron
(veja `casos.txt` e o relatório em `../relatorio/`):

| Fase | Casos | Instantes do FS | N |
|---|---|---|---|
| 1 (6.1) | referência, I = I_ref (padrão 31,05 mm/h, onde o GIM dá Γ = 1,333, o do Ceron) | 0, 6, 12, 18 e 24 h | 5 000 |
| 2 (6.2) | CoV(k_s) 0; 0,75; 0,9; 1 · CoV(c') 0,1; 0,5; 0,7 · CoV(φ') 0,05; 0,15; 0,2 | 24 h (k) ou 0 e 24 h | 2 000 |
| 2 (6.3) | s = 1,5; 2; 5; 10; 20; 400 (L_x = 20s m, L_y = 2s m) | 0 e 24 h | 2 000 |
| 2 (6.4) | α = k_h/k_v = 2, 3, 4, 5 | 24 h | 2 000 |
| 3 (6.5) | I = 24, 27, 29, 33, 36 e 45 mm/h | 24 h | 2 000 |

**Custo estimado** (realização em serial na estação: 75 s só com 24 h, 117 s com 0 e 24 h, 240 s com os cinco
instantes; núcleo do fermi ~2,6 vezes mais lento):

| Fase | Horas de CPU no fermi | Com 4 nós × 12 processos |
|---|---|---|
| 1 | ~870 h | ~18 h |
| 2 | ~2 900 h | ~60 h |
| 3 | ~650 h | ~14 h |
| Total | ~4 400 h | ~3,8 dias |

Com `PARALELO=24` (hyperthreading) deve render mais por nó; não foi medido.

## Rodar

Do PC, com a VPN ligada (cada comando pede a senha do fermi).

### 0. Mandar o código

```bash
cd ~/projects/neopz-master && git push origin UnsaturatedSoil
```

### 1. Código no fermi

Se `~/neopz-master` ainda existe lá:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil'
```

Se foi apagado:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'git clone -b UnsaturatedSoil https://github.com/diogocecilio/neopz-master.git neopz-master'
```

### 2. Compilar no nó de login

A primeira vez leva de 15 a 30 min. O comando roda em segundo plano, então pode desconectar:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && nohup ./compilar.sh > compilar.log 2>&1 &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 3 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/compilar.log'
```

A última linha deve ser `Compilado: .../bin/SlopeRainMonteCarlo (<commit>_<data>)`.

### 3. Teste rápido

Um job de ~15 min: o determinístico e 2 realizações da referência.

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && sbatch teste.sbatch'
ssh diogo@fermi.cesup.ufrgs.br 'cat neopz-master/Projects/SlopeRandomFieldRefinement/cluster/logs/chuvaTeste_*.out'
```

Com I = 31,05 mm/h, o determinístico deve dar FS(0 h) ≈ 2,094 e FS(24 h) ≈ 1,126.

### 4. Submeter

Use 4 nós com 12 processos cada (o limite é de 4 jobs rodando por usuário):

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && ./submeter.sh 4'
```

Para começar só pela referência, que dá ~18 h:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CASOS=fase1 ./submeter.sh 4'
```

### 5. Acompanhar

```bash
ssh diogo@fermi.cesup.ufrgs.br 'neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
```

Mostra os jobs, as últimas linhas dos logs e, por caso, as tarefas feitas, as realizações e o P_f contado em 24 h.

### 6. Baixar e gerar o relatório, no PC

```bash
~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

- Pede a senha uma vez e baixa `$DATA/chuva_mc` para `~/projects/runs/fermi_chuva`.
- Roda `../posproc/dados_chuva.py` e `relatorio_chuva.py`.
- O PDF vai para `relatorio_parcial/` enquanto há jobs na fila e para `relatorio/` depois.

### Parar, retomar, estender

- **Parar** sem perder nada (os processos terminam a tarefa em curso):

  ```bash
  ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/chuva_mc/PARAR'
  ```

  Para cortar na hora: `scancel -n chuvaMC`. As tarefas em curso ficam travadas e são refeitas na próxima submissão
  sem jobs na fila.
- **Retomar:** `./submeter.sh 4` de novo. As tarefas feitas são puladas, e as interrompidas são liberadas e refeitas.
  Também dá para submeter mais jobs com outros ainda rodando: todos dividem a mesma fila.
- **Estender:** `ULTIMA=9999 CASOS=ref ./submeter.sh 4` leva a referência a 10 000 realizações. O Ceron usou blocos de
  10 000 até CoV(P_f) < 5 %. Com P_f ≈ 1 %, isso pede ~40 000 realizações; com 10 %, ~3 600.

## Como funciona

- **Tarefa:** um bloco de `BLOCO` realizações (padrão 20, umas 2 h no fermi) de um caso. Cada uma é uma rodada do
  `SlopeRainMonteCarlo` em `$SAIDA/<caso>`.
- **Ordem:** o bloco 0 de todos os casos, depois o 1, ..., para que os casos andem juntos.
- **Divisão entre processos:** cada processo (`tarefas.sh`; o `job.sbatch` roda `PARALELO` por nó) pega a próxima
  tarefa livre criando o arquivo `$SAIDA/<caso>/.tarefas/<a>_<z>` com o `noclobber` do bash (`O_EXCL`), que só um
  processo consegue, também entre nós. Não é `mkdir`: o coreutils em Rust (uutils) de Ubuntu recente não cria
  diretórios de forma atômica, e no teste dois processos pegavam a mesma tarefa.
- **Fim de uma tarefa:** a tarefa está feita quando existe `mc_mapa_<a>_<z>.csv`, gravado só no fim da rodada.
- **Bloco 0 de cada caso:** roda também o determinístico (`det=1`) e grava os campos dele (`mc_real-1/`), de onde sai a
  figura da banda do mecanismo.
- **`$SAIDA/campanha.txt`:** guarda `IREF` e `BLOCO`, que definem os arquivos. `submeter.sh` recusa outros valores
  na mesma `SAIDA`.
- **Variáveis de `submeter.sh`:** `SAIDA` (padrão `$DATA/chuva_mc`), `CASOS` (`todos`, ids ou fases), `IREF`,
  `BLOCO`, `ULTIMA`, `PARALELO`, `TEMPO` (padrão 2 dias) e `PARTICAO`.

## Arquivos de cada caso (`$SAIDA/<caso>/`)

| Arquivo | Conteúdo |
|---|---|
| `mc_<a>_<z>.csv` | uma linha por realização e instante (colunas abaixo) |
| `mc_mapa_<a>_<z>.csv` | quantas realizações têm a banda do mecanismo em cada célula de 1 m, e quantas dessas falham |
| `mc_<a>_<z>.log` | a saída do programa |
| `mc_parametros.txt` | solo, chuva, campos, malhas e instantes do caso |
| `mc_real-1/` | o determinístico: arquivos do `LearningMohrCoulomb2/animacao.py`, `campos.csv` e `mecanismo.csv` |

Colunas de `mc_<a>_<z>.csv`; as que existem no `mc_resumo.csv` do rebaixamento têm o mesmo nome:

| Coluna | Conteúdo |
|---|---|
| `real`, `t_h` | realização (−1: determinístico) e instante (h) |
| `FS_min`, `situacao` | FS do SRM no instante; `ok`, `transiente` (o transiente parou antes) ou `fator` |
| `niveis` | o FS de cada nível do refinamento, separados por `;` |
| `banda_*` | a banda de ruptura (√J2(ε_p) ≥ 0,1 do máximo): área, profundidade, extensão, centroide, p_w, c', φ' e k médios |
| `massa_*` | a massa que desliza (\|Δu\| desde o estado geostático ≥ 0,5 do máximo): área, centroide, c', φ', k |
| `c_med`, `phi_med`, `k_mult_med` | médias dos campos no solo todo |
| `nivel_x30`, `nivel_x35`, `nivel_x40` | cota do nível d'água na crista, no meio da face e no pé (categorias I a IV do Rahardjo) |
| `pw_min`, `pw_max`, `nos_face` | faixa de p_w e nós da superfície na face de percolação |
| `eq_*`, `s_*` | equações e tempos do transiente e do fator |
| parâmetros | `solo`, `I`, `chuva`, `hw`, `metodo`, `cov_c`, `cov_phi`, `cov_k`, `Lx`, `Ly`, `alfa`, `KL_variancia` |

As variáveis do Spearman da Seção 6.1 são `banda_c`, `banda_phi_deg` e `massa_k_mult`, como no rebaixamento.

## Teste local

Rodado antes da submissão, numa pasta de teste:
- 2 processos, blocos de 1 realização, 2 casos;
- uma tarefa interrompida e retomada.

Resultados:
- os processos pegam tarefas diferentes;
- a tarefa interrompida fica travada até a limpeza do `submeter.sh` e depois é refeita;
- o determinístico sai só no bloco 0;
- o `status.sh` mostra a tabela.
