# Rodar os contraexemplos A3–A6 no seu computador

Divisão do trabalho: o ambiente remoto roda o A2 (fase de ar, `scripts/fase_ar.py`) e o A1 (dilatância,
`SlopeRainfallPV`); o seu computador roda os itens **A3, A4, A5 e A6** (99 casos independentes, modelo do artigo:
percolação desacoplada e p_a = 0). Cada caso usa 1 núcleo e leva 10–20 min; com N núcleos o total é ≈ 99 × 15 / N min
(≈ 1,5 h com 16 núcleos).

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
pip install numpy matplotlib
mkdir runs && cd runs
python <neopz>/Projects/SlopeRainfall/scripts/contraexemplos.py rodar --itens A3 A4 A5 A6 \
       --jobs <número de núcleos físicos> --exe <build>/Projects/SlopeRainfall/SlopeRainfall
```
Os resultados vão para `runs/resultados/contraexemplos/<caso>/`. Casos com `summary.txt` são pulados: se a execução
for interrompida, basta rodar o mesmo comando de novo. Para ver o progresso:
`ls resultados/contraexemplos/*/summary.txt | wc -l` (total 99).

As figuras preliminares podem ser geradas no próprio computador:
```
python <neopz>/Projects/SlopeRainfall/scripts/contraexemplos.py figuras
```

## 4. Devolver os resultados

Compacte só os arquivos pequenos (os `snaps_*.csv` têm ≈ 0,8 MB por caso e não são necessários para as figuras):
```
cd runs
zip -r contraexemplos.zip resultados/contraexemplos -x "*/snaps_*" "*/log.txt"
```
e anexe o `contraexemplos.zip` na conversa (ou faça commit de `resultados/contraexemplos` em
`Projects/SlopeRainfall/resultados/` no branch). Os resultados remotos (A1, A2) e os seus serão reunidos no relatório.
