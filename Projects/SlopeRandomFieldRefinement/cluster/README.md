# Monte Carlo no fermi: chuva (Rahardjo) e rebaixamento do Ceron

Esta pasta roda duas campanhas, escolhidas pela variável `CAMPANHA` (veja `campanha.sh`):
- `chuva` (padrão): o Monte Carlo da chuva, descrito logo abaixo;
- `ceron`: o rebaixamento rápido do Ceron com o transiente, com Monte Carlo direto e amostragem por importância
  (seção [Rebaixamento do Ceron com amostragem por importância](#rebaixamento-do-ceron-com-amostragem-por-importância));
- `chuva_ai`: a chuva com amostragem por importância, com o FS de 6 em 6 h até 48 h
  (seção [Chuva com amostragem por importância](#chuva-com-amostragem-por-importância));
- `solos`: os três solos do Rahardjo com I = 1 k_s, Monte Carlo direto com o FS de 6 em 6 h até 48 h
  (seção [Os três solos](#os-três-solos)).

As duas usam a mesma fila de tarefas (`tarefas.sh`), a mesma compilação e os mesmos comandos, com
`CAMPANHA=ceron` na frente para a segunda.

## Chuva

O `montecarlo.cpp` (alvo `SlopeRainMonteCarlo`) roda o Monte Carlo do talude de Rahardjo et al. (2007) sob chuva:
- talude da Fig. 2 do Rahardjo, com H = 10 m a 45°, na `TriGMesh`;
- solo f50: k_s = 1e-5 m/s, Fredlund-Xing com a = 50 kPa, c' = 10 kPa, φ' = φ^b = 26°;
- transiente da chuva com esqueleto rígido e ar passivo, como no SEEP/W do artigo;
- FS pelo SRM, com refinamento adaptativo em h refeito do zero em cada instante (`TriGMesh(1)` + 4 níveis);
- c', φ' e k_s lognormais, gerados por Karhunen-Loève com M = 2000 e as estatísticas da Tabela 2 de Vargas Ceron et al.
  (2025).

A realização n tem os mesmos números aleatórios em todos os casos. Os casos são o análogo das Seções 6.1 a 6.5 do Ceron
(veja `casos.txt` e o relatório em `Projects/PlanoDeTrabalho/campanhas/`):

| Fase | Casos | Instantes do FS | N |
|---|---|---|---|
| 1 (6.1) | referência, I = I_ref (padrão 31,05 mm/h, onde o GIM dá Γ = 1,3345, perto do 1,336 do Ceron) | 0, 6, 12, 18 e 24 h | 5 000 |
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

As horas com 4 nós são a CPU dividida por 48. No relógio, a fase 1 leva ~21 h: são 250 tarefas de ~3,5 h para 48
processos. Com `PARALELO=24` (hyperthreading) deve render mais por nó; não foi medido.

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

Para começar só pela referência, que dá ~21 h:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CASOS=fase1 ./submeter.sh 4'
```

### 5. Acompanhar

```bash
ssh diogo@fermi.cesup.ufrgs.br 'neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
```

Mostra os jobs, as últimas linhas dos logs e, por caso, as tarefas feitas, as realizações e o P_f do menor FS dos
instantes (em 24 h nos casos com 0 e 24 h).

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

## Rebaixamento do Ceron com amostragem por importância

O `rebaixamento_mc.cpp` (alvo `SlopeDrawdownMonteCarlo`) roda o Monte Carlo do rebaixamento rápido de Vargas Ceron et
al. (2025), Seção 6.1:
- talude de 5 m a 45°;
- transiente do `main.cpp`: o reservatório desce da crista ao pé em 0,1 h, com k_v = 1e-6 m/s;
- Γ pelo GIM no fim do rebaixamento, com refinamento adaptativo (`TriGMesh(1)` + 3 níveis);
- c', φ' e k_v lognormais da Tabela 2 do Ceron (KL com M = 2000).

O determinístico dá Γ = 1,3253 (o artigo dá 1,336). No PC, com 486 realizações, P_f ≈ 18,5% (Ceron: 11,5%).

Casos (`casos_ceron.txt`):

| Caso | N | O que é |
|---|---|---|
| `ref_direto` | 10 000 | Monte Carlo direto: o piloto das amostragens e a referência para conferir que elas não têm viés |
| `ref_ai_banda` | 5 000 | amostragem por importância com a média dos ξ deslocada na direção da banda de ruptura do determinístico (um escalar por campo) |
| `ref_ai_modos` | 5 000 | amostragem por importância com o deslocamento nos 20 primeiros modos de cada campo (`UnsaturatedSoil/formulacao/form_MIS`) |

Como funciona a amostragem por importância:
- **Piloto:** os dois casos `ref_ai_*` só começam quando as realizações 0 a 999 do `ref_direto` acabaram. Até lá, os
  processos adiantam o direto.
- **Deslocamento:** cada bloco calcula a média dos ξ das realizações do piloto que falharam (entropia cruzada), sorteando
  de novo os ξ delas pela semente. O resultado é o mesmo em todos os blocos e fica em `mc_desloc_<a>_<z>.txt`.
- **Peso:** cada realização grava o peso exato na coluna `log_peso`. P_f = Σ peso × falha / N não tem viés.
- **Independência:** as sementes são diferentes, então as três amostras são independentes.

**Custo:** ~30 s por realização no PC, ~80 s no fermi. As 20 000 realizações são ~440 h de CPU, ou ~10 h com 4 nós.

**Comandos**, do PC, com a VPN ligada. Os passos 0 a 2 são os da chuva e compilam os dois programas:

```bash
cd ~/projects/neopz-master && git push origin UnsaturatedSoil
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil'
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && nohup ./compilar.sh > compilar.log 2>&1 &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 3 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/compilar.log'
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CAMPANHA=ceron sbatch teste.sbatch'
ssh diogo@fermi.cesup.ufrgs.br 'cat neopz-master/Projects/SlopeRandomFieldRefinement/cluster/logs/chuvaTeste_*.out'
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CAMPANHA=ceron ./submeter.sh 4'
ssh diogo@fermi.cesup.ufrgs.br 'CAMPANHA=ceron neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
CAMPANHA=ceron ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

- **Compilar:** a última linha do log deve ser `Compilado: .../bin/SlopeRainMonteCarlo e SlopeDrawdownMonteCarlo (...)`.
- **Teste:** ~5 min; o determinístico deve dar Γ ≈ 1,325.
- **Status:** P_f ponderado de cada caso e o CoV dele.
- **Baixar:** baixa `$DATA/ceron_ai` para `~/projects/runs/fermi_ceron_ai` e grava em `resultado.txt` o P_f de cada
  amostragem contra o direto, o CoV e quantas realizações diretas dariam o mesmo CoV.

Parar e retomar como na chuva: `touch /mnt/data/diogo/ceron_ai/PARAR`, `scancel -n ceronMC`, ou
`CAMPANHA=ceron ./submeter.sh 4` de novo.

**Mais casos (Seções 6.2 a 6.5):** acrescente em `casos_ceron.txt` um par direto + amostragem para cada caso. O direto
também é o piloto, por exemplo:

```
covc=0.1_direto  2  2000  covc=0.1
covc=0.1_ai      2  2000  covc=0.1 piloto=@SAIDA@/covc=0.1_direto piloto_n=1000 ai=banda banda=@SAIDA@/covc=0.1_direto/mc_real-1/mecanismo.csv semente=1
```

Depois, submeta de novo: a fila pega os casos novos.

Com P_f pequeno, o piloto tem poucas falhas. Nesse caso, o conjunto elite são os 10% com menor Γ (opção `elite`): é um
nível intermediário da entropia cruzada, com um deslocamento menor que o ótimo, mas ainda sem viés.

## Os três solos

O `casos_solos.txt` roda os três solos da Tabela 1 do Rahardjo com a chuva da Fig. 5: I = 1 k_s por 24 h (f10 a
360 mm/h, f50 a 36 mm/h, f100 a 3,6 mm/h).
- **Amostra:** Monte Carlo direto, com 1000 realizações por solo e o FS de 6 em 6 h até 48 h. O FS em 0 h fica só no
  f50, porque no determinístico ele é o mesmo nos três.
- **Comparação pareada:** a realização n tem os mesmos campos nos três solos.
- **Determinístico:** o f10 fica abaixo de 1 de ~6 h a ~26 h, o mínimo do f50 é 1,035 em 24 h e o do f100 é 1,886 em
  32 h. Por isso o P_f vai de quase 100% (f10) a quase 0 (f100). O que se compara é o FS(t) com a sua dispersão.

**Custo:** ~22 min por realização no fermi; as 3000 realizações levam ~1 dia com 4 nós. Os blocos têm 5 realizações.

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CAMPANHA=solos ./submeter.sh 4'
ssh diogo@fermi.cesup.ufrgs.br 'CAMPANHA=solos neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
CAMPANHA=solos ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

O `baixar.sh` grava em `~/projects/runs/fermi_chuva_solos/resultado.txt` o resumo de cada solo
(`posproc/resumo_tempos.py`):
- em cada instante: o FS médio, o CoV, P(FS < 1) e a P_f do ajuste lognormal;
- no menor FS dos instantes: a P_f e o instante em que o mínimo ocorre.

## Chuva com amostragem por importância

A referência da chuva (f50, I = 31,05 mm/h) tem P_f ≈ 1% em 24 h, pelos pilotos. Com isso, o Monte Carlo direto vê
poucas falhas: 5000 realizações dão CoV(P_f) ≈ 14%. A amostragem por importância (`montecarlo.cpp` com `piloto=`; a
lógica está em `ImportanceSampling.h`) precisa de bem menos realizações (Seção 7 do `Projects/PlanoDeTrabalho/campanhas/relatorio_chuva.pdf`).

Por que até 48 h:
- dependendo do k, o FS mínimo não ocorre em 24 h: com k baixo, ele vem depois da chuva (no f100 determinístico, às
  32 h);
- a falha passa a ser FS < 1 em algum dos instantes 6, 12, ..., 48 h (o mínimo no tempo);
- o peso de cada realização não depende do instante: a mesma amostra dá a P_f(t) e a P_f do mínimo.

Casos (`casos_chuva_ai.txt`):

| Caso | N | O que é |
|---|---|---|
| `ref48_piloto` | 200 | Monte Carlo direto, o piloto. Com poucas falhas, o conjunto elite são os 10% de menor FS |
| `ref48_ai` | 200 | amostragem por importância com o deslocamento na direção da banda do determinístico em 24 h; começa quando o piloto acaba |

**Custo:** ~17 min por realização no fermi (o transiente até 48 h e 8 instantes do SRM). Os blocos têm 1 realização,
para os 48 processos trabalharem; as 400 realizações levam ~2,5 h com 4 nós.

**Precisão esperada** com P_f ≈ 1%: CoV(P_f) de ~14% (estado limite linear) a ~25% com 200 realizações. Para estender:
`ULTIMA=999 CASOS=ref48_ai CAMPANHA=chuva_ai ./submeter.sh 4` leva a amostragem a 1000 realizações.

**Comandos:** compile como na chuva e:

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && CAMPANHA=chuva_ai ./submeter.sh 4'
ssh diogo@fermi.cesup.ufrgs.br 'CAMPANHA=chuva_ai neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
CAMPANHA=chuva_ai ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

O `status.sh` mostra a P_f do mínimo no tempo, já com os pesos. O `baixar.sh` grava em
`~/projects/runs/fermi_chuva_ai/resultado.txt` a P_f do mínimo e a de cada instante, contra o piloto.

**Com a campanha do Ceron rodando:** o fermi deixa 4 jobs rodando por usuário, então os jobs `chuvaAI` esperam na fila.
Para recompilar só o `SlopeRainMonteCarlo` com jobs `ceronMC` na fila, use `FORCAR=1 ./compilar.sh`. Isso é seguro
enquanto a NeoPZ não muda: a `libpz.so` dos jobs fica igual, e os executáveis são trocados de forma atômica.

## Etapas 0 e 1 do plano de trabalho (`CAMPANHA=etapa1`)

O plano está em `Projects/PlanoDeTrabalho/README.md`. Os casos ficam em `casos_etapa1.txt`, todos com o SRM novo
(`srm_f0`, padrão 0,8: o caminho de F começa abaixo de 1):
- **Etapa 0:** o Cho (2010) que rodava no PC:
  - Monte Carlo direto com r = −0,5 (4000) e r = −0,7 (1000);
  - a AI pela banda (1000 cada);
  - a AI com o ponto de projeto do FORM de 10 termos por campo (`r05_ai_form`, 1000).
- **Etapa 1:** o SRM antigo (`srm_f0=1`) contra o novo nas mesmas realizações:
  - as 1000 primeiras do Cho2010 (`r05_antigo`);
  - o f50 da chuva em 24 h, 600 em cada versão.

O primeiro job roda antes da fila os FORM com M = 10 e 20 (`form.sh`, com o `Projects/Cho2010/form_kl.py`; 12
avaliações ao mesmo tempo, ~1,5 h) e depois entra na fila. O `r05_ai_form` espera o ponto de projeto
(`espera=@SAIDA@/form_M10/ustar.txt`). Ao todo são ~290 h de núcleo, umas 6 a 7 h com 4 nós.

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "python3 --version; ./compilar.sh && CAMPANHA=etapa1 ./submeter.sh 4" > etapa1.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 8 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/etapa1.log'
ssh diogo@fermi.cesup.ufrgs.br 'CAMPANHA=etapa1 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
CAMPANHA=etapa1 ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

O primeiro comando compila e submete em segundo plano: pode desconectar. O `status.sh` mostra também a última iteração
de cada FORM. `FORM_M=` (vazio) no `submeter.sh` não roda os FORM, e `FORM_M="10"` roda só um. Sem `python3` no nó,
o FORM grava o `.ERRO`, o `r05_ai_form` fica de fora e o resto roda.

## Etapa 3b do plano de trabalho (`CAMPANHA=etapa3b`)

Quantas realizações a ISMC (amostragem por importância no ponto de projeto do FORM) precisa para CoV(P_f) = 5% com o
campo de Nyström de 2000 modos (`casos_etapa3b.txt`; o plano em `Projects/PlanoDeTrabalho/README.md`):

- **Os 6 FORM** rodam antes da fila em **todos** os jobs (`FORM_LISTA` do `campanha.sh`), e as travas os dividem:
  - `form_r05_M100` e `form_r07_M100`: diferenças centrais com M = 100, os mais caros, primeiro;
  - `form_r05_M20` e `form_r07_M20`: diferenças centrais com M = 20;
  - `form_r05_gim20` e `form_r07_gim20`: GIM com o gradiente modal.
  Eles leem a base com os 100 primeiros modos, `$SAIDA/pontos/base_nys2000_M100.txt`.
- **As 30 ISMC:** 3 centros × 2 valores de r × 5 sementes, com os 2000 modos (`$SAIDA/pontos/base_nys2000.txt`) e
  `cov_alvo=0.05`. Cada uma espera o seu FORM (`espera=`).
- **A conferência de versão** (`conf_r07_s3`): as 10 primeiras realizações do `r07_ai_form20` da etapa3.

A base (525 MB em texto, mais a cópia de 27 MB) vai do PC para o fermi uma vez só, antes de submeter. O rsync
comprime no caminho:

```bash
rsync -avz --progress --rsync-path="mkdir -p /mnt/data/diogo/etapa3b/pontos && rsync" ~/projects/runs/cho2010/nystrom2000/base_nys2000.txt ~/projects/runs/cho2010/nystrom2000/base_nys2000_M100.txt diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/etapa3b/pontos/
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar.sh && CAMPANHA=etapa3b ./submeter.sh 4" > etapa3b.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 8 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/etapa3b.log; CAMPANHA=etapa3b neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh'
CAMPANHA=etapa3b ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

O `baixar.sh` roda o `Projects/Cho2010/ismc_parada.py`. Os diretos de referência são os já baixados das etapas 1 e 3:
baixe o `r07_direto` de 110 000 antes, com `CAMPANHA=etapa3 baixar.sh`.

## Tabelas 4 e 5 com o GIM, pareado com o SRM (`CAMPANHA=etapa3b_gim` e `ls_gim.sbatch`)

Os casos do Cho (2010) das Tabelas 4 e 5 do `relatorio_comparacao.pdf` (r = -0,5 e -0,7) com o fator pelo GIM
(`metodo=gim`, Γ; falha: Γ < 1) em vez do SRM (FS). Base da comparação análise limite (`la_cho.sbatch`) × MEF SRM ×
MEF GIM.

- **Pareado:** as opções e as sementes são as das rodadas SRM, só com `metodo=gim` a mais. A realização n depende só
  da semente, de n, de M e da base, então é o mesmo campo (mesmos `c_med`, `phi_med` e `log_peso` no CSV) e Γ e FS
  se comparam realização a realização. Uma semente só por caso.
- **Sem FORM novo:** a ISMC usa os centros u* da etapa 3b (`$DATA/etapa3b/form_*/ustar.txt`) e a base de Nyström de
  lá (`$DATA/etapa3b/pontos/base_nys2000.txt`). O line sampling usa as direções do SRM: `form_r0x_M20` (já no fermi) e
  `form_r0x_gim2000` (do PC: suba os dois `form.json`, primeiro comando abaixo).
- **Sem recompilar:** o `bin/Cho2010` do fermi (commit ac5aabf) já tem `metodo=gim`, e o código do Cho2010 não mudou
  depois.

Casos (`casos_etapa3b_gim.txt`, blocos de 10 como nas referências):

| Caso | O que é | Referência SRM | Análises |
|---|---|---|---|
| `conf_srm_r05` | SRM, realizações 0 a 9 (conferência: FS iguais) | `etapa1/r05_direto` | 10 |
| `conf_srm_r07` | SRM, realizações 0 a 9 (conferência: FS iguais) | `etapa3/r07_direto` | 10 |
| `r0x_nys_<M20,M100,gim20>_s21` | ISMC GIM, semente 21, CoV 5%, no mínimo 2000 (`cov_n=2000`) | `etapa3b/<mesmo id>` | ~12 300 |
| `ls_gim.sbatch` | LS GIM, semente 7, direções M20 e gim2000, CoV 5%, 100 a 600 linhas | `~/projects/runs/etapa3b_ls/ls` | ~1 800 |

O `cov_n=2000` faz a ISMC GIM rodar pelo menos as 2000 realizações da ISMC da análise limite e as k da parada do SRM
(1270 a 2120). O N em que o CoV cruza 5% sai depois, pelo `ismc_parada.py cov_n=500` (a regra do SRM), que o
`baixar.sh` roda com os FORM da etapa 3b (custo com o FORM) e o z contra o direto GIM. O `NMIN=100` do LS cobre as
100 linhas do SRM e da análise limite. O `ls_gim.sbatch` copia a base de Nyström (525 MB) para o disco local do nó,
porque cada análise do LS é um processo novo que a lê.

Sem Monte Carlo direto com o GIM (decisão do usuário, 2026-11-20): a referência da P_f é o direto SRM, nas mesmas
realizações. Custo (o SRM leva ~66 a 69 s por análise no fermi; o GIM, pelo PC, não mais): ~14 100 análises, ~270 h
de núcleo. Com 2 jobs da campanha e o `lsGim` (3 rodando; o 3º job da campanha espera o `lsGim` acabar, com
`DEPENDE`), dá umas 5 a 7 h. Fica 1 vaga rodando para o `laCho`.

Antes de submeter, confira no fermi a versão do `bin/Cho2010` (tem de ser a ac5aabf) e o `random` do python3 (o Z do
LS tem de ser o mesmo do PC: `-0.2558802884476004 0.511431512516514`):

```bash
cd ~/projects/runs/etapa3b_ls && rsync -avR --rsync-path="mkdir -p /mnt/data/diogo/etapa3b_gim && rsync" form_r05_gim2000/form.json form_r07_gim2000/form.json diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/etapa3b_gim/
ssh diogo@fermi.cesup.ufrgs.br 'cat neopz-master/Projects/SlopeRandomFieldRefinement/cluster/bin/versao.txt; python3 -I -c "import random;r=random.Random(7);print(repr(r.gauss(0.,1.)),repr(r.gauss(0.,1.)))"'
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && CAMPANHA=etapa3b_gim ./submeter.sh 2 && J=$(sbatch --parsable ls_gim.sbatch) && J=${J%%;*} && echo "lsGim $J" && CAMPANHA=etapa3b_gim DEPENDE=$J ./submeter.sh 1'
ssh diogo@fermi.cesup.ufrgs.br 'squeue -u diogo; CAMPANHA=etapa3b_gim neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh; tail -n 4 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/logs/lsGim_*.out'
CAMPANHA=etapa3b_gim ~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar.sh
```

O `baixar.sh` baixa para `~/projects/runs/fermi_etapa3b_gim` (com o `ls/`) e faz o `resultado.txt`: a conferência,
cada caso GIM contra a sua referência SRM (`Projects/Cho2010/pareado_gim.py mc`: campo igual?, falhas e P_f de cada
um, pares discordantes, k de Γ ~ FS^k, e a tabela `pareado.csv` com real, FS_srm, Gamma_gim), o LS linha a linha
(`pareado_gim.py ls`) e o `ismc_parada.py`. As referências SRM têm de estar baixadas em `~/projects/runs`.

- **Γ = 0** (GIM sem equilíbrio nem com 1% da gravidade, `situacao` `fator`) conta como falha no `tarefas.sh`, no
  `status.sh`, no `amostragem_importancia.py`, no `compara_pares.py` e no `line_sampling.py`.

## Tabelas 4 e 5 com a análise limite do Ceron (`la_cho.sbatch`)

Os mesmos casos do Cho (r = -0,5 e -0,7) com a análise limite cinemática do Ceron (2025) em Python
(`Projects/CeronCinematico/cho`), nos **mesmos campos** do MEF: o gerador, a KL analítica e a base de Nyström são
cópias do C++, então a realização n da LA é a realização n do MEF (conferido: log_peso igual, campos até 5e-6).
Grade dos campos na LA: h = 0,025 m (com 0,1 o Γ sai ~2,5 % baixo). Uma semente por caso (ISMC 21, LS 7).

| Etapa | O que roda | Custo (estimado) |
|---|---|---|
| `checa` | python (venv), random.gauss, determinístico (Γ 1,777, FS 1,203), caminhos do MEF no `$DATA`, Ψ de Nyström (cache 3,9 GB, uma vez) e o tempo de 48 realizações | ~0,5 h |
| `direto R="r05 r00 r07"` | MC direto LA até CoV 5 % (r05 ≥ 24 000, r07 ≥ 110 000: todas as do MEF, pareadas; r00 = r 0, só para o Ceron, 5,50 %). Γ em todas (~1,6 s); FS da LA (~13 s) nas 2400 primeiras e perto da falha | ~150 h de CPU |
| `amostragem R="r05 r07"` | FORM LA (M 20 e 100), LS LA (100 linhas, semente 7) nas direções do MEF (M20, gim2000) e da LA, ISMC LA (semente 21, 2500) nos centros do MEF (M20, M100, gim20) e da LA | ~60 h de CPU |
| `pareado` | (opcional) LA e `bin/Cho2010` nos mesmos xi, para conferir o caminho xi= | ~5 h |

No fermi (uma vez: `$HOME/venv_la/bin/pip install numpy scipy matplotlib`):

```bash
cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster
J=$(ETAPA=checa sbatch --parsable la_cho.sbatch)
ETAPA=direto R="r05 r00 r07" sbatch --dependency=afterok:$J la_cho.sbatch
ETAPA=amostragem R="r05 r07" sbatch --dependency=afterok:$J la_cho.sbatch
```

Saídas em `$DATA/cho_la/h0.025` (direto/, form/, ls/, ismc/); o fim de cada job imprime a tabela LA × MEF
(`resumo.py`). Baixar: `rsync -av --exclude cache diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/cho_la/ ~/projects/runs/fermi_cho_la/`.

## Line sampling do Ceron com a agua do artigo (`contorno=ceron`)

O `ls_ceron.sbatch` calcula a P_f do rebaixamento do Ceron (2025) com os contornos do artigo (`contorno=ceron
agua=sat trebaix=0.0001`: percolação saturada depois do rebaixamento instantâneo; o artigo dá P_f = 11,50%):
- FORM com o gradiente modal e M = 2000 termos por campo (`form_kl.py prog=rebaixamento`, 1 análise por iteração, até 15);
- line sampling na direção do FORM, com K = 2000, parando quando CoV(P_f) ≤ 5% (mínimo 30, máximo 600 linhas, 24 em
  voo);
- Γ no permanente: `tfim=20 dt=5` (4 passos; o mesmo Γ que `dt=0.1` até com k = 0,05 k_v no solo todo);
- malhas: transiente `tri3`, fator `tri1` com 5 níveis (com 3 níveis o Γ fica ~4% alto e a P_f baixa).

Usa uma build e um executável próprios (`build_ceron/` e `bin_ceron/`, do `compilar_ceron.sh`). Assim pode compilar e
rodar com jobs das campanhas na fila. Resultados em `$DATA/ceron_ls` (`form/form.json`, `ls/ls.json`, `ls/linhas.csv`, `ls/perdidas.csv`).
Linhas perdidas (transiente sem convergir ou Γ NaN) ficam fora da conta, vão para o `perdidas.csv` e são somadas no
`ls.json` (`fracao_perdida`); o resumo avisa quando passam de 1% (a P_f pode sair baixa). Γ = 0 (GIM sem equilíbrio)
continua como falha.
Tempo estimado: 0,5 a 1,5 h para o FORM e 1 a 3 h para as linhas (4 a 5 h se chegar a 600).

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar_ceron.sh && sbatch ls_ceron.sbatch" > ceron_ls.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master/Projects/SlopeRandomFieldRefinement/cluster && tail -n 3 ceron_ls.log; squeue -u diogo; tail -n 6 /mnt/data/diogo/ceron_ls/form.log; tail -n 4 logs/lsCeron_*.out'
rsync -av diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/ceron_ls/ ~/projects/runs/ceron_ls/
```

- **Parar** sem perder linhas: `ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/ceron_ls/ls/PARAR'`.
- **Retomar:** `sbatch ls_ceron.sbatch` de novo (apague antes o `PARAR`, se houver). O FORM pronto é pulado e as
  linhas feitas não são refeitas.
- **Outro alvo:** `COV=0.03 N=1500 sbatch ls_ceron.sbatch` continua as mesmas linhas até o novo CoV.

## Como funciona

- **Tarefa:** um bloco de `BLOCO` realizações (padrão 20, umas 2 h no fermi) de um caso. Cada uma é uma rodada do
  `SlopeRainMonteCarlo` em `$SAIDA/<caso>`.
- **Ordem:** o bloco 0 de todos os casos, depois o 1, ..., para que os casos andem juntos.
- **Divisão entre processos:** cada processo (`tarefas.sh`; o `job.sbatch` roda `PARALELO` por nó) pega a próxima
  tarefa livre criando o arquivo `$SAIDA/<caso>/.tarefas/<a>_<z>` com o `noclobber` do bash (`O_EXCL`), que só um
  processo consegue, também entre nós. Não é `mkdir`: o coreutils em Rust (uutils) de Ubuntu recente não cria
  diretórios de forma atômica, e no teste dois processos pegavam a mesma tarefa.
- **Fim de uma tarefa:** a tarefa está feita quando existe `mc_mapa_<a>_<z>.csv`, gravado só no fim da rodada. Se a
  rodada termina com erro, a trava ganha uma linha `ERRO`. A tarefa fica parada até a limpeza do `submeter.sh`, sem
  jobs na fila, e então é refeita.
- **Dependência:** um caso com `piloto=@SAIDA@/<id> piloto_n=<n>` só começa quando os blocos das realizações 0..n−1 do
  caso `<id>` acabaram. Um caso com `espera=<arquivo>,<arquivo>` só começa quando os arquivos existem. Se um deles tem
  um `<arquivo>.ERRO` (o `form.sh` grava um quando o FORM falha), o caso fica de fora. Sem outra tarefa livre, o
  processo espera 60 s por vez, até `ESPERA_MAX` (3 h).
- **Opções da fila**, que não vão para o programa:
  - `prog=<nome>` troca o programa do caso por outro de `bin/`;
  - `espera=` (acima);
  - `ref=<id>` e `par=<id>` marcam o direto de referência e a comparação pareada para o pós-processamento.
  - `cov_alvo=<cv>` para o caso quando o CoV da P_f ponderada chega a `cv`. A conta usa as realizações já prontas,
    em ordem, com pelo menos `cov_n` (500) realizações e `cov_falhas` (100) falhas sorteadas. O caso grava
    `ATINGIU` e não lança mais blocos; o N do caso vira só o teto.
- **P_f no `status.sh`:** cada realização conta com o menor FS dos seus instantes, e só se todos tiverem fator. Na chuva
  com 0 e 24 h, é o FS de 24 h.
- **Opções:** nas opções dos casos, `@SAIDA@` vira a pasta da campanha e `IREF` vira a intensidade de referência.
- **Bloco 0 de cada caso:** roda também o determinístico (`det=1`) e grava os campos dele (`mc_real-1/`), de onde sai a
  figura da banda do mecanismo.
- **`$SAIDA/campanha.txt`:** guarda `IREF` e `BLOCO`, que definem os arquivos. `submeter.sh` recusa outros valores
  na mesma `SAIDA`.
- **Variáveis de `submeter.sh`:** `SAIDA` (padrão `$DATA/chuva_mc`), `CASOS` (`todos`, ids ou fases), `IREF`,
  `BLOCO`, `ULTIMA`, `PARALELO`, `TEMPO` (padrão 2 dias), `PARTICAO` e `DEPENDE` (o id de um job: os novos só começam
  quando ele acabar).

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
