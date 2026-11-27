# Rahardjo2007

Reprodução da **Série B** de H. Rahardjo, T. H. Ong, R. B. Rezaur e E. C. Leong (2007), *Factors controlling
instability of homogeneous soil slopes under rainfall*, J. Geotech. Geoenviron. Eng. 133(12):1532–1543.

- Talude homogêneo de 10 m a 45° (Fig. 2), nível d'água inicial 5 m abaixo do pé subindo 7° para dentro do maciço,
  sucção inicial limitada a 75 kPa.
- Três solos (Tabela 1): f10,-4, f50,-5 e f100,-6 (Fredlund-Xing com a = 10, 50 e 100 kPa, n = m = 1, θs = 0,45;
  k_w = k_s Θ⁴; c' = 10 kPa, φ' = φ^b = 26°, γ = 20 kN/m³).
- Sete intensidades por solo, chuva de 24 h, FS até 48 h (Fig. 5); FS mínimo contra a intensidade (Figs. 6 e 7, com mais
  7 pontos que só aparecem nelas); a posição do nível d'água no FS mínimo (Tabelas 3 e 5).
- Modelos do transiente (`modelo=`):
  - **(iii)**, o padrão: esqueleto rígido e ar passivo, só a percolação da água (o SEEP/W do artigo);
  - **(ii)**: acoplado (esqueleto elástico, E = 20 MPa e ν = 0,3, que o artigo não dá), ar passivo;
  - **(i)**: acoplado e com o ar ativo (p_a = 0 só na superfície; k_ra de Brooks–Corey, λ = 3, piso 1e-4).
- FS pelo SRM (o artigo usa o Bishop simplificado do SLOPE/W), com as pressões do transiente congeladas.
- Resultado do (iii): reproduz o artigo dentro do erro de leitura das figuras (relatório em `relatorio/`). O (i) é a
  base para uma Eq. 6 corrigida para o FS com ar ativo.

## Executável

`Rahardjo2007 [solo=f10|f50|f100] [I=<mm/h>] [modelo=iii|ii|i] [E=2e4] [nu=0.3] [tempos=artigo|<k>h|t1,t2,...]
[horas=48] [dt=0.1] [nref=4] [dx=1] [eta0=0.003] [r=1.15] [etamax=0.025] [saida=rah] [so_transiente=0|1]
[grava_malha=<prefixo>] [campos=t1,t2,...]` (veja o cabeçalho do `main.cpp`). Saídas: `<pre>_fs.csv` e
`<pre>_nivel.csv`, com `<pre>` = `<saida>_<solo>_<I>` no (iii) e `<saida>_<modelo>_<solo>_<I>` nos outros.

## Campanha no fermi (`cluster/`)

`cluster/casos.txt` tem as 28 rodadas do (iii), e `cluster/casos_i.txt` e `cluster/casos_ii.txt` as mesmas com o (i) e
o (ii) (FS de hora em hora), 7 por nó e 3 threads cada (`rah.sbatch`): os 4 nós rodam todas de uma vez. A conferência de malha
(`cluster/casos_malha.txt`) e os campos do (iii) (`cluster/casos_campos.txt`) rodam no PC com o mesmo `rah.sbatch`
(`SAIDA=... EXE=... CASOS_ARQ=... PARALELO=3 THREADS=4 ./rah.sbatch`).

No fermi a campanha usa uma cópia própria do código (`~/neopz-rah`, um worktree do `~/neopz-master`), para não mexer
nos arquivos das outras campanhas. O projeto vai do PC para `~/rah_src` por rsync e é copiado para lá antes de compilar.
Um comando por vez, no terminal do PC:

```bash
# 1. mandar o projeto
rsync -a ~/projects/neopz-master/Projects/Rahardjo2007 ~/projects/neopz-master/Projects/CMakeLists.txt diogo@fermi.cesup.ufrgs.br:rah_src/
# 2. copiar, compilar (so o que mudou) e submeter o (i) em $DATA/rahardjo2007_i
ssh diogo@fermi.cesup.ufrgs.br 'nohup setsid bash -c "set -x; cp -a ~/rah_src/Rahardjo2007 ~/rah_src/CMakeLists.txt ~/neopz-rah/Projects/ && cd ~/neopz-rah/Projects/Rahardjo2007/cluster && JOBS=12 ./compilar.sh && SAIDA=/mnt/data/diogo/rahardjo2007_i CASOS_ARQ=\$PWD/casos_i.txt ./submeter.sh 4" > ~/rah_lanca_i.log 2>&1 < /dev/null & echo lancado'
# 3. acompanhar
ssh diogo@fermi.cesup.ufrgs.br 'tail -n 3 ~/rah_lanca_i.log; squeue -u diogo; cd ~/neopz-rah/Projects/Rahardjo2007/cluster && SAIDA=/mnt/data/diogo/rahardjo2007_i CASOS_ARQ=casos_i.txt ./rah.sbatch resumo | tail -n 8'
# 4. o (ii) em $DATA/rahardjo2007_ii, com o mesmo executavel (o casos_ii.txt e o casos_i.txt com modelo=ii). Foi
#    lancado junto com a volta do f100_I5.4 do (i), que ficou com a trava vazia (veja abaixo): o primeiro submeter.sh,
#    sem job na fila, libera a trava
ssh diogo@fermi.cesup.ufrgs.br 'cd ~/neopz-rah/Projects/Rahardjo2007/cluster && sed "s/modelo=i\b/modelo=ii/" casos_i.txt > casos_ii.txt && SAIDA=/mnt/data/diogo/rahardjo2007_i CASOS_ARQ=$PWD/casos_i.txt PARALELO=1 THREADS=12 ./submeter.sh 1 && SAIDA=/mnt/data/diogo/rahardjo2007_ii CASOS_ARQ=$PWD/casos_ii.txt ./submeter.sh 4 && squeue -u diogo'
```

Trava órfã: no (i), o `.trava` do f100_I5.4 foi criado vazio no instante em que os 28 processos pegavam casos (um soluço do
sistema de arquivos do fermi), e nenhum processo o rodou. O `resumo` agora mostra esses casos como "trava vazia"; o
`submeter.sh`, sem job `rah2007` na fila, libera as travas dos casos sem `fim`.

Download no PC: `cluster/baixar.sh` (o iii, para `~/projects/runs/rahardjo2007/fermi`), `cluster/baixar.sh i` (para
`fermi_i`) e `cluster/baixar.sh ii` (para `fermi_ii`). Relatório: `relatorio/gera.sh`.
