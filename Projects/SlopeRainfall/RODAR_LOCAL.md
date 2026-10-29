# Rodar os contraexemplos no seu computador

Os itens **A3, A4, A5 e A6** (99 casos independentes, modelo do artigo: percolação desacoplada e p_a = 0) já rodaram
(resultados em `resultados/contraexemplos/`). Faltam o **A1** (dilatância, 15 casos, `SlopeRainfallPV`, modelo
acoplado com ar ativo) e o **A2** (fase de ar, 22 casos, `scripts/fase_ar.py`): veja a seção 5. Cada caso usa 1
núcleo; os do A3–A6 levam 10–40 min, os acoplados do A1 e do A2, 30–150 min.

## 1. Código

```
git fetch origin
git checkout claude/relaxed-volta-1r4gog
git pull
```

## 2. Compilar (Linux, macOS ou WSL)

```
mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON \
      -DUSING_MKL=ON -DMKL_ROOT=<pasta do MKL> <caminho do neopz>
ninja SlopeRainfall
```
* Com MKL o sistema linear usa Pardiso (bem mais rápido). Sem MKL omita `-DUSING_MKL` e `-DMKL_ROOT` (usa a skyline,
  mais lenta). O MKL pode vir do oneAPI ou do `pip install mkl-devel` (neste caso `MKL_ROOT` é o prefixo do Python,
  p. ex. `/usr/local`, e talvez seja preciso criar os links `libmkl_*.so` → `libmkl_*.so.2`).
* Teste rápido: `./Projects/SlopeRainfall/SlopeRainfall ttotal=1 tfs=1 out=teste` (≈ 2 min).

## 3. Rodar

```
pip install numpy matplotlib      # (use python3 se o comando python não existir)
mkdir runs && cd runs
python3 <neopz>/Projects/SlopeRainfall/scripts/contraexemplos.py rodar --itens A3 A4 A5 A6 \
       --jobs <número de núcleos físicos> --exe <build>/Projects/SlopeRainfall/SlopeRainfall
```
Os resultados vão para `runs/resultados/contraexemplos/<caso>/`. Casos com `summary.txt` são pulados: se a execução
for interrompida, basta rodar o mesmo comando de novo. Para ver o progresso:
`ls resultados/contraexemplos/*/summary.txt | wc -l` (total 99).

As figuras preliminares podem ser geradas no próprio computador:
```
python3 <neopz>/Projects/SlopeRainfall/scripts/contraexemplos.py figuras
```

## 4. Devolver os resultados

Compacte só os arquivos pequenos (os `snaps_*.csv` têm ≈ 0,8 MB por caso e não são necessários para as figuras):
```
cd runs
zip -r contraexemplos.zip resultados/contraexemplos -x "*/snaps_*" "*/log.txt"
```
e anexe o `contraexemplos.zip` na conversa (ou faça commit de `resultados/contraexemplos` em
`Projects/SlopeRainfall/resultados/` no branch).

## 5. A1 e A2

`scripts/rodar_A1_A2.sh` compila `SlopeRainfall` e `SlopeRainfallPV` e lança os dois lotes em segundo plano (um caso
por núcleo; continuam com o terminal fechado). Na pasta dos resultados:
```
cd ~/projects/runs
<neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh rodar <build> 4 4     # núcleos do A1 e do A2
<neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh status                # casos prontos (de 15 e 22) e em execução
<neopz>/Projects/SlopeRainfall/scripts/rodar_A1_A2.sh zip                   # no fim: a1_a2.zip para anexar
```
Logs em `lote_A1.log` e `lote_A2.log`. Um lote que já está rodando não é lançado de novo e os casos com
`summary.txt` são pulados: se a máquina for desligada, basta repetir o `rodar`. O total é de ≈ 4–5 h com 8 núcleos
(o A1 é o mais longo).

## 6. Benchmarks (Zhang, Tang, Pradhan) e patamar inferior do FS

`scripts/rodar_benchmarks.py` roda os 29 casos (zhang 8, patamar 6, tang 13, pradhan 2), os mais longos primeiro, um
núcleo por caso, e pula os que já têm `summary.txt` (depois de uma interrupção basta repetir o `rodar`). Na pasta dos
resultados:
```
cd ~/projects/runs
cmake --build <build> --target SlopeRainfall -j 8
S=<neopz>/Projects/SlopeRainfall/scripts
nohup python3 $S/rodar_benchmarks.py rodar --exe <build>/Projects/SlopeRainfall/SlopeRainfall --jobs 8 > benchmarks.log 2>&1 &
python3 $S/rodar_benchmarks.py status                                  # prontos e em execução, por grupo
python3 $S/benchmarks.py --runs benchmarks --figuras figuras           # figuras e tabelas
zip -qr benchmarks.zip benchmarks                                      # no fim: para anexar
```
Para rodar só alguns grupos: `... rodar --exe ... pradhan patamar`; `lista` mostra os comandos de cada caso. Os mais
longos são os dois casos de Pradhan (271 h de chuva em malhas de ≈ 6200 triângulos) e o patamar `f100_iso`
(6144 h de superfície alagada).
