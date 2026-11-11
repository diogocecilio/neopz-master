#!/bin/bash
# compilar.sh: compila a NeoPZ e o alvo SlopeMonteCarloCluster em ./build e copia o executavel para ./bin, com a versao
# (o commit do git, "+mod" se havia alteracoes locais, e a data) em bin/versao.txt, que vai para a coluna versao de cada
# realizacao. No fermi, no no de login (os nos de computacao nao tem cmake), com o MKL carregado:
#     source ../../LearningMohrCoulomb2/cluster/modulos.sh && ./compilar.sh
# A primeira compilacao (a biblioteca inteira) leva uns 15 a 30 min; depois, so o que mudou.
# Ambiente: USAR_MKL (ON: PARDISO), MKL_THREAD_MODEL (seq), JOBS (4), BUILD (./build) e FORCAR=1 (compila mesmo com jobs
# slopeMC na fila: eles usam a libpz.so de build/, que seria trocada no meio da conta).
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
BUILD="${BUILD:-$AQUI/build}"
USAR_MKL="${USAR_MKL:-ON}"
if [ -z "${FORCAR:-}" ] && command -v squeue > /dev/null; then
    fila="$(squeue -h -u "${USER:-$(id -un)}" -n slopeMC -o %i)" || { echo "o squeue nao respondeu: tente de novo (ou FORCAR=1)" >&2; exit 1; }
    [ -z "$fila" ] || { echo "ha jobs slopeMC na fila usando a libpz.so de $BUILD: espere terminarem (ou FORCAR=1)" >&2; exit 1; }
fi
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL="$USAR_MKL")
if [ "$USAR_MKL" = ON ]; then ARGS+=(-DMKL_THREAD_MODEL="${MKL_THREAD_MODEL:-seq}"); else ARGS+=(-DUSING_LAPACK=OFF); fi
# Catch2 local, se o pacote cluster de LearningMohrCoulomb2 ja o baixou (configurar sem internet)
CATCH2="$AQUI/../../LearningMohrCoulomb2/cluster/deps/Catch2"
[ -f "$CATCH2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$CATCH2")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target SlopeMonteCarloCluster -j "${JOBS:-4}"
mkdir -p "$AQUI/bin"
cp "$BUILD/Projects/SlopeMonteCarlo/SlopeMonteCarloCluster" "$AQUI/bin/.SlopeMonteCarloCluster.novo"
mv -f "$AQUI/bin/.SlopeMonteCarloCluster.novo" "$AQUI/bin/SlopeMonteCarloCluster" # troca atomica
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2>/dev/null || echo sem_git)"
[ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=no 2>/dev/null)" ] && v="$v+mod"
echo "${v}_$(date '+%F_%T')" > "$AQUI/bin/versao.txt"
echo "Compilado: $AQUI/bin/SlopeMonteCarloCluster ($(cat "$AQUI/bin/versao.txt"))"
