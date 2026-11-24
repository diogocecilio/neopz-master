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
- **Baixar:** baixa `$DATA/ceron_ai` para `~/projects/runs/_arquivo/2026-11-14_fermi_ceron_ai/fermi_ceron_ai` e grava em `resultado.txt` o P_f de cada
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

## Ceron 6.1 a 6.5 com a agua do artigo: FORM modal + line sampling (`ls_ceron_casos.sbatch`)

As Secoes 6.1 a 6.5 de Vargas Ceron et al. (2025) com o codigo corrigido (a agua do artigo, `contorno=ceron agua=sat
trebaix=0.0001`, como no `ls_ceron.sbatch`), os 26 casos das Tabelas 3 a 6 (`casos_ls_ceron.txt`):

| Secao | Casos |
|---|---|
| 6.1 | referencia |
| 6.2 | CoV(k_v) 0; 0,75; 0,9; 1 · CoV(c') 0,1; 0,5; 0,7 · CoV(phi') 0,05; 0,15; 0,2 |
| 6.3 | s = 1,5; 2; 5; 10; 20; 400 (L_x = 20 s m, L_y = 2 s m) |
| 6.4 | alfa = k_h/k_v = 2, 3, 4, 5 |
| 6.5 | h_w/H = 0,5; 0,6; 0,7; 0,8; 0,9 |

Em cada caso, em tres fases (os jobs dividem tudo por travas em `$SAIDA`):
1. **FORM** com o gradiente modal, M = 2000 termos por campo (1 analise por iteracao, ate 15) e o **deterministico**
   (as medias, `det=1`): cada job roda ate 24 dessas tarefas ao mesmo tempo;
2. **line sampling** na direcao do FORM (K = 2000, semente 7 em todos os casos), ate CoV(P_f) <= 5% (30 a 600
   linhas, 24 em voo), um caso por vez;
3. **Monte Carlo direto pareado**, 500 realizacoes por caso (a realizacao n tem os mesmos campos em todos os casos,
   como no artigo), so para as densidades, a media e o CoV de Gamma das figuras; roda por ultimo (ou enquanto um
   job espera FORM de outro).

Base comum: `tfim=20 dt=5` (Gamma no permanente), transiente `tri3`, fator `tri1` com 5 niveis. Os deterministicos
no PC (21/11) ja batem com o artigo: ref 1,351 (1,336), alfa 2 a 5: 1,552 / 1,692 / 1,799 / 1,885 (1,533 / 1,674 /
1,783 / 1,872), h_w/H 0,5 a 0,9: 1,687 / 1,505 / 1,390 / 1,328 / 1,312 (1,671 / 1,494 / 1,383 / 1,322 / 1,307).

**Custo estimado** (a referencia no PC: FORM em 4 analises, LS de 153 linhas e 363 analises para CoV 5%; ~2,5 min por
analise no fermi com 24 processos por no): ~6 analises de FORM, ~360 a 1400 de LS e 500 de MC por caso, ~30 000
analises ao todo, umas 12 a 16 h com 4 nos (o LS de todos os casos em ~8 h; o MC depois). Menos MC: `NMC=300`.

```bash
cd ~/projects/neopz-master && git push origin UnsaturatedSoil
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar_ceron.sh && ./submeter_ls_ceron.sh 4" > ls_ceron_casos.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 5 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/ls_ceron_casos.log; neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status_ls_ceron.sh'
~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar_ls_ceron.sh
```

- O `compilar_ceron.sh` recompila so o que mudou (o C++ nao mudou desde o `ls_ceron.sbatch`) e recusa compilar com
  jobs `lsCeron` ou `lsCeronC` na fila. O `submeter_ls_ceron.sh` confere que o `bin_ceron/` tem o `contorno=ceron`.
- `status_ls_ceron.sh`: os jobs, os ultimos eventos (`$SAIDA/eventos.log`) e a tabela por caso (deterministico,
  FORM, LS parcial, MC) contra o artigo.
- `baixar_ls_ceron.sh` (PC): baixa para `~/projects/runs/fermi_ceron_ls_casos`, grava `resultado.txt` e
  `resultado.csv` e faz as figuras da Secao 6 do artigo em `figs/` (`figuras_ls_ceron.py`: Figs. 25, 27 a 33 e o
  deterministico, com o Ceron tracejado; P_f do line sampling com IC 95%, do FORM e do MC; densidades, media e CoV
  de Gamma do MC).
- Parar sem perder nada: `ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/ceron_ls_casos/PARAR'` (as linhas e
  os blocos em curso terminam). Voltar: apagar o `PARAR` e `./submeter_ls_ceron.sh` de novo (as travas dos jobs que
  morreram sao liberadas; o LS continua do `linhas.csv`). Com outras campanhas rodando, submeta menos nos
  (`./submeter_ls_ceron.sh 2`) ou com `DEPENDE=<job>`.

## Caso A: Hu et al. (2024), FORM modal + line sampling + Monte Carlo (`hu2024_casos.sbatch`)

O caso principal da verificacao da chuva (relatorio `Exemplos probabilisticos de chuva nao saturada`, Secao "A) Caso
principal: Hu et al. (2024)"): talude de 10 m a 45 graus num dominio de 70 x 20 m, nivel d'agua horizontal na cota 5 m
com a succao inicial limitada a 50 kPa, lateral esquerda e base fechadas, lateral direita com carga de 5 m abaixo do
nivel, chuva por metro inclinado em toda a superficie (com a face de percolacao onde empocaria), van Genuchten-Mualem da
Fig. 4 (theta_s 0,399, theta_r 0,024, alfa 0,0357 1/kPa, n 1,71), chi = S_e, gamma 17,8 kN/m3, c' lognormal 8 kPa (CoV
0,15), phi' lognormal 24 graus (CoV 0,10), r(c, phi) = -0,5, k_s = 2e-5 m/s fixo. Os campos sao de Nystrom no proprio
talude (nucleo exp(-|dx|/20 - |dy|/2), grade 0,3 x 0,1 m, 34 950 nos, M = 150 modos por campo, 86,9% da variancia). O FS
e o SRM com as pressoes congeladas (fator `tri1` + 4 niveis); o transiente e o da `tri4`; a falha e o menor FS abaixo
de 1 (primeira passagem).

Os 6 cenarios de `casos_hu2024.txt` (os 5 do relatorio e o I30_D24, em que o MC direto confere o LS):

| caso | chuva | P_f do artigo | FS deterministico no fim da chuva (tri3) | N do MC | teto do LS (linhas) |
|---|---|---|---|---|---|
| I30_D12 | 30 mm/h por 12 h | 0,012 | 1,296 | 500 | 60 |
| I20_D24 | 20 mm/h por 24 h | 0,030 | 1,226 | 500 | 60 |
| I40_D12 | 40 mm/h por 12 h | 0,091 | 1,132 | 4000 | 600 |
| I20_D36 | 20 mm/h por 36 h | 0,117 | 1,118 | 4000 | 600 |
| I50_D24 | 50 mm/h por 24 h | 0,7781 | 0,908 | 300 | 600 |
| I30_D24 | 30 mm/h por 24 h | 0,221 | 1,065 | 2000 | 600 |

O N do MC nao e o do CoV de 10% com a P_f do artigo: no piloto de 22/11 (10 realizacoes, 50 mm/h por 24 h) o desvio do
FS e ~0,37 do de Hu (0,063, 0,034 e 0,047 em 0, 24 e 72 h, contra 0,165, 0,091 e 0,130; as medias batem: 1,524, 0,905 e
1,178 contra 1,58, 0,93 e 1,25), entao as nossas P_f ficam longe das do artigo. O FORM do I30_D12 no PC deu beta ~ 12
(nao convergiu em 10 iteracoes: o FS e muito nao linear tao longe da origem). Veja o cabecalho do `casos_hu2024.txt`.

Com k_s fixo, o transiente e o mesmo em todas as realizacoes: cada caso o calcula uma vez (`so_transiente=1
grava_pressoes=`) e todas as analises o leem (`pressoes=`). Fases (os jobs dividem tudo por travas em `$SAIDA`):
0. a base de Nystrom (com o python do venv, `$HOME/venv_la/bin/python`; uns 2 min) e o transiente de cada caso (`tri4`,
   uns 25 min no fermi com 4 threads);
1. o deterministico (FS de 6 em 6 h de 0 a 72 h) e o FORM com o gradiente modal (M = 150, ate 10 iteracoes de 1 analise);
2. o line sampling na direcao do FORM ate CoV(P_f) <= 5% (30 a 600 linhas, 24 em voo), um caso por vez;
3. o Monte Carlo direto (as N realizacoes do caso, em blocos de 10, intercalados) e a conferencia `mc13/`: as 20 primeiras
   realizacoes de cada caso com o FS nos 13 instantes, que mede em quantas o menor FS cai no fim da chuva.

Instantes (`TEMPOS=fim`, o padrao): o FS no fim da chuva D, com refinamento adaptativo, e em D + 6 h, na malha de D. O
deterministico de todos os cenarios e as 10 realizacoes do piloto dao o menor FS em t = D (o artigo: "the FOS reaches
minimum value when the rainfall stop"); `TEMPOS=todos` faz os 13 instantes em toda analise (13 SRM em vez de ~1,3).

**Custo estimado:** uma analise (o FS em D e D + 6) leva ~24 s em serial no PC (~60 a 95 s por processo no fermi, com 24
por no); ~11 300 do MC, ~1 200 equivalentes da conferencia e ~1 500 do LS: ~14 000 analises, umas 4 a 5 h com 4 nos
(mais ~30 min do transiente e da base no inicio).

```bash
cd ~/projects/neopz-master && git push origin UnsaturatedSoil
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar_hu2024.sh && ./submeter_hu2024.sh 4" > hu2024.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 5 neopz-master/Projects/SlopeRandomFieldRefinement/cluster/hu2024.log; neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status_hu2024.sh'
~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar_hu2024.sh
```

- `KS=5.2e-6` (o patamar da Fig. 4b) roda a mesma campanha com o outro k_s, numa pasta propria (`$DATA/hu2024_ks5.2e-6`):
  `KS=5.2e-6 ./submeter_hu2024.sh 4` (e `KS=5.2e-6` nos outros comandos). O deterministico de 22/11 favorece o 2e-5
  (FS em 24 h a 50 mm/h: 0,908 contra 1,104 com o 5,2e-6; o artigo tem 0,93).
- `status_hu2024.sh`: os jobs, os ultimos eventos e a tabela por caso (deterministico, FORM, LS parcial, MC, conferencia)
  contra o artigo; `baixar_hu2024.sh` (PC) baixa para `~/projects/runs/fermi_hu2024_ks2e-5` e grava `resultado.txt`,
  `resultado.csv` e `resumo.json`.
- Parar sem perder nada: `ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/hu2024_ks2e-5/PARAR'`. Voltar: apagar o
  `PARAR` e `./submeter_hu2024.sh` de novo. Com outras campanhas rodando, submeta menos nos ou com `DEPENDE=<job>`.

## Caso B: Liu et al. (2024), permanente com campos gaussianos (`liu2024_casos.sbatch`)

O passo intermediario da verificacao (relatorio `Exemplos probabilisticos de chuva nao saturada`, Secao "B) Passo
intermediario"): o talude de Griffiths e Fenton (5 m a 1:2; pontos (0;0), (30;0), (30;5), (20;5), (10;10), (0;10)),
cargas de 6 m a esquerda e 1 m a direita abaixo do nivel, base e laterais acima do nivel fechadas, percolacao permanente
com q = 1e-7 m/s (0,36 mm/h) na superficie (por metro inclinado, o fluxo unitario do GeoStudio; a Fig. 1 desenha setas
verticais, e a chuva vertical muda o FS em +0,1%), Fredlund-Xing (a 5 kPa, n 2, m 1, theta_s 0,4) com o k pela integral de
Fredlund (Eq. 9, com C(psi) e psi_r = 3000 kPa; psi_aev = 2,6 kPa, a tangente), k_s 1e-5 m/s (CoV 0,6), c' 10 kPa (CoV
0,3), phi' 15 graus (CoV 0,2), r(c, phi) = -0,5, phi^b = 2 phi'/3 por ponto (chi = tan(2phi'/3)/tan(phi')), gamma 20.
Campos de Nystrom com o nucleo gaussiano exp(-dx^2/30^2 - dy^2/3^2) (alfa = 0), grade 0,2 x 0,1 m (11 225 nos), M = 100
modos por campo (variancia retida 1,0000: o gaussiano cai rapido; lambda_k da grade de 0,1 x 0,05 m dentro de 4% ate o
modo 100, em que lambda_100/lambda_1 = 7e-14). O permanente e a malha `tri3` (0,625 m); o fator, `tri1` + 3 niveis.

Codigo novo (`montecarlo.cpp`, `RainTransient`, `TPZUnsaturatedHydraulics.h`; os padroes antigos ficam bit a bit):
- `q0=` (mm/h): o estado inicial e a percolacao permanente com esse fluxo (`RainTransient::SteadyState`: Euler implicito
  com dt dobrando de 0,1 h ate 1e8 h; para quando a maior mudanca de p_w num passo com dt >= 1e4 h fica abaixo de 1e-6 kPa;
  no caso B, 19 passos e 9 s);
- `fx=a,n,m,theta_s` e `kfredlund=<psi_aev>` (`psi_r=`): o k_r de Fredlund tabelado contra S_w (o novo
  `TPZUnsaturatedHydraulics::ETabulated`, opcional; o material `TPZMatUnsaturatedUPwPa` nao muda);
- `chi=ponto` com `phib=` ou `phib_razao=`: chi(x) = tan(phi^b)/tan(phi'(x)) nas pressoes congeladas do fator.

Verificacao no PC (22/11, `~/projects/runs/liu2024`):
- nivel d'agua do permanente: y = 5,14 m na secao X-X (x = 10) e 4,56 m na Y-Y (x = 15), contra ~5,0 e ~4,5 na Fig. 3;
  p_w na base 49,9 e 44,1 kPa (Fig. 3: 50 e ~44); succao no topo 13,2 kPa contra ~15 (o patamar de k_r = q/k_s);
- FS deterministico 1,5764 (`tri3`, 4 niveis; 1,5785 com 3 niveis, 1,5783 na `tri4`, 1,5758 com o fator na `tri2`) contra
  1,525 do Bishop do artigo (+3,4%); sem a parcela da succao (phi^b = 0), 1,518;
- FORM modal do cenario 2 (4 iteracoes no PC): beta ~ 4,75 a 4,88 (o FS oscila ~0,003 entre iteracoes), contra 3,95 da
  Tabela 4; phi' pesa 0,78 e c' 0,22;
- piloto do cenario 3 (4 realizacoes): FS 1,40 a 1,51; uma analise leva 1,5 a 4 min em serial no PC (o SRM deste talude
  e lento: 5 000 a 16 000 equacoes, caminho longo).

Os 2 cenarios de `casos_liu2024.txt` (alfa = 0, q/k_s = 0,01):

| cenario | campos | P_f do artigo (Tabela 4) | N do MC | teto do LS (linhas) |
|---|---|---|---|---|
| S2 | c', phi' (k_s fixo: um permanente so, lido por todas as analises) | 4,08e-5 (beta 3,95) | 2 500 | 600 |
| S3 | k_s, c', phi' (um permanente por analise) | 1,19e-4 (beta 3,68) | 10 000 | 600 |

O MC do S3 confere a verificacao da Tabela 2 (4,0e-3 por hipercubo latino com 2500 amostras, a mesma de Jiang et al.
2014): com 10 000 realizacoes, ~40 falhas se a P_f for 4e-3 e ~1 se for 1,19e-4. As Tabelas 2 e 4 do artigo nao fecham
(~40 vezes) para "o mesmo" modelo; o resumo mostra as duas, sem forcar nenhuma. No S3 o gradiente modal da 0 nos termos de
k_s (o fator congela as pressoes): a direcao do FORM fica em c' e phi', e o LS continua sem vies (o artigo da beta = 15 so
com k_s aleatorio). Fases: 0. a base (gauss, ~20 s) e o permanente do S2; 1. o deterministico e o FORM (M = 100, ate 10
iteracoes); 2. o LS ate CoV 5% (30 a 600 linhas); 3. o MC em blocos de 10.

**Custo estimado:** uma analise leva 1,5 a 4 min em serial no PC (~4 a 10 min por processo no fermi, com 24 por no): MC
~12 500 analises, LS ~2 x 600 x 2,5, FORM ~20: ~16 000 analises, ~1,6 a 4 mil horas de CPU, ou ~17 a 40 h com 4 nos
(NMC_MAX menor encurta; o LS do S2 sai primeiro).

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar_liu.sh && ./submeter_liu2024.sh 2" > liu2024.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status_liu2024.sh'
~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar_liu2024.sh
```

- Com os jobs do caso A (huA) ou do Ceron na fila: `./submeter_liu2024.sh 2` (2 nos) ou `DEPENDE=<id do ultimo huA>
  ./submeter_liu2024.sh 4` (os 4 esperam aquele job terminar). O fermi roda 4 jobs por usuario e aceita 8 na fila.
- Parar sem perder nada: `ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/liu2024/PARAR'`.

## Caso C: Liu e Wang (2023), chuva aleatoria e CRPC (`liuwang2023_casos.sbatch`)

O terceiro caso (relatorio, Secao "C) Terceiro caso"): talude de 10 m num dominio de 74 x 30 m (pontos (0;0), (74;0),
(74;20), (43;20), (31;30), (0;30)), cargas de 20 m e 10 m abaixo do nivel inicial reto de (0;20) a (74;10), base e
laterais acima do nivel fechadas (suposicao: o artigo nao diz), van Genuchten-Mualem (theta_s 0,358, theta_r 0, a = 6,711
kPa: alfa 0,149 1/kPa, n 1,289), k_s 3e-6 m/s, phi^b = 15 graus fixo (chi = tan15/tan(phi'(x)) por ponto), gamma 20, c'
5 kPa, phi' 32 graus, campos lognormais sem correlacao cruzada, nucleo exponencial 20/2 m por Nystrom (grade 0,5 x 0,1 m,
37 000 nos: lambda_1..200 dentro de 0,56% dos da grade de 0,2 x 0,1 m; M = 200, 84,7% da variancia). Estado inicial: o
permanente com q0 = chuva anual / 8760 h, a chuva anual normal (1200 mm, CoV 0,1). Chuva por metro inclinado (o fluxo
unitario do GeoStudio, setas normais a face na Fig. 5).

**Ainda nao submeter (22/11): o transiente nao converge neste solo.** Com van Genuchten n = 1,289 o Newton do transiente
entra num ciclo de 2 iteracoes (|R_p| alternando ~1e-4, a busca linear nao acha descida): no solo das medias ele para
entre ~27 h (I = 4 mm/h) e ~70 h (I = 36 mm/h), quando a superficie toda vira face de percolacao, e com os campos ja em 2
a 3 h de chuva (piloto do MC no PC: 1 de 2 realizacoes parou). Com n = 2 o mesmo transiente passa de 60 h; com n = 1,5
para em 45 h. Nao resolveram: passo subdividido ate 0,00625 h (`dt_min=`), tolerancia maior e mais iteracoes (`newton=`),
busca linear com 12 metades, k_r de Mualem suavizado perto da saturacao e uma retencao com curvatura finita na saturacao
(testes em `~/projects/runs/liuwang2023/det`; as duas ultimas nao ficaram no codigo). O `submeter_liuwang2023.sh` so
submete com `FORCAR=1`.

**A chuva da Tabela 2 nao fecha com a P_F1 do artigo.** O artigo imprime sigma_D = -0,30 h (impossivel; Jun et al. 2017
nao esta acessivel); com |sigma_D| = 0,30 h, P(D > 24 h) = 4,3e-4, P(D > 48 h) = 6,7e-5 e P(D > 65 h) = 3,1e-5. As CRPC do
artigo tem D_C >= 65 h (o parametro c vai de 65,4 a 188,8 h), entao a P_F1 dessa chuva ficaria abaixo de ~3e-5, contra
4,6e-4 a 4,8e-3 impressos. No nosso solo das medias o FS0 e 1,271 e a 36 mm/h o FS ainda esta >= 1 em 48 h, entao a
nossa P_F1 tambem fica abaixo de ~7e-5 com essa chuva. O `gpd=` troca os parametros (ex.: o sigma_D certo, quando houver).

Codigo e campanha (prontos, com a ressalva acima):
- `chuva_aleatoria=1` (I, D da Pareto generalizada com copula de Frank, q0 normal; FS em 0 e em D; `<saida>_chuva_*.csv`);
- `crpc=100,36,15,9,6,4,3,2,1` (CRPC por realizacao: FS em 6, 12, 24, 48, 96, 150 e 200 h ate o primeiro < 1, bissecao
  nos instantes de 1 em 1 h, D_C por interpolacao; as I menores param depois de uma que nao derruba ate 200 h);
- `crpc_liuwang.py`: o ajuste D_C = a I^b + c, I_lb, P_1 pela integral da densidade de (I, D) (1D na condicional da
  copula, Simpson) e P_F1 = media de P_1 com o CoV da Eq. 14; e o MC direto com o descarte de FS0 < 1;
- `liuwang2023_casos.sbatch` (jobs liuC): o solo fixo (CRPC, uma tarefa por intensidade) e o MC direto dos niveis baixo e
  intermediario (5500 cada, o N do artigo), em blocos de 10; `submeter/status/baixar_liuwang2023.sh`,
  `resumo_liuwang2023.py`. Uma analise do MC leva ~2 a 4 min no PC (dois SRM de ~40 s + o permanente + o transiente).

## Cho (2014): coluna 1D sob chuva, Monte Carlo + FORM/LS (`cho2014_casos.sbatch`)

A coluna de talude infinito de Cho (2014, Eng. Geol. 171: 11-20) no `SlopeRainColumn` (`coluna.cpp`; spec dos passos 1
a 3, Secao D): H = 2 m, k_s lognormal 1D (base de Nystrom do segmento, exp(-|dz|/lv), M = 100), van Genuchten-Mualem
com o k_r tabelado (eps 1e-3), I = 10,8 mm/h por 24 h sobre succao inicial de 20 kPa, FS de Lu e Godt depois de cada
passo (beta 40 graus, c' 5 kPa, phi' 32 graus). Os padroes de producao do executavel (estudos V2.x em
`~/projects/runs/coluna_estudos`): quad com massa consistente, ny 80, dt 0,1 h, newton=1e-6,80,1, skyline.

- `casos_cho2014.txt`: os 15 casos da Fig. 10 (CoV 0,5/1,0/1,5 x lv 0,5/2/4/8/16 m) com N = 5000 e as P_f do artigo
  em 18 e 24 h; os estudos V2.5 (ny 200), V2.7 (eps 1e-4) e V3.3 (tri concentrada e consistente), pareados; as linhas
  Ippisch (Q1) comentadas.
- `casos_cho2014_ls.txt`: F1a/F1b (c05_l0.5 em 18 h, falha fs/tempo), F2 (16 h, fs) e F3 (c05_l16 em 20 h, tempo):
  `form_kl.py prog=coluna gradiente=dif h=0.5 algoritmo=ihlrf` (M = 100), `line_sampling.py` ate CoV 5% e
  `bissecao_ls.py`, depois do MC do caso (a referencia).
- Fases do job: deterministicos; MC em blocos de 100 intercalados entre os casos; FORM/LS/bissecao de cada linha do LS
  com o MC do caso pronto. As bases (5, uma por lv) sao feitas pelo `submeter_cho2014.sh` no no de login (segundos).
- `resumo_cho2014.py` (status e baixar): P_f no instante e de primeira passagem com Wilson e os limites da censura, z
  contra a Fig. 10 e contra a referencia em Python (`~/projects/runs/cho2014_python`, no PC), os alvos das Figs. 3, 5, 6,
  7 e 9, e o LS contra o MC do mesmo caso (z, cos(alfa_FORM, xi medio das falhas), alfa^2, custo).

**Custo estimado:** ~7 s por realizacao no PC (quad ny 80), ~2,6 vezes no fermi: ~530 h de CPU no MC (5,5 h com 4 nos,
11 h com 2), mais ~0,5 h por FORM (15 iteracoes de 201 analises em 24 processos).

```bash
ssh diogo@fermi.cesup.ufrgs.br 'cd neopz-master && git pull origin UnsaturatedSoil && cd Projects/SlopeRandomFieldRefinement/cluster && nohup setsid bash -c "./compilar_coluna.sh && ./submeter_cho2014.sh 4" > cho2014.log 2>&1 < /dev/null &'
ssh diogo@fermi.cesup.ufrgs.br 'neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status_cho2014.sh'
~/projects/neopz-master/Projects/SlopeRandomFieldRefinement/cluster/baixar_cho2014.sh
```

- Com os jobs liuB rodando, os jobs cho14 alem dos 4 por usuario esperam na fila (o fermi aceita 8).
- Parar sem perder nada: `ssh diogo@fermi.cesup.ufrgs.br 'touch /mnt/data/diogo/cho2014/PARAR'`.
