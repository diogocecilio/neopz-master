#!/bin/bash
# compilar.sh: compila a NeoPZ e o SlopeRainMonteCarlo numa pasta de build PROPRIA (padrao: ./build), sem tocar na build/
# do SlopeMonteCarlo, e copia o executavel para ./bin com a versao (o commit do git, "+mod" se havia alteracoes locais, e
# a data) em bin/versao.txt. Os jobs usam bin/ (a troca e atomica), mas a libpz.so continua em build/: nao compile com
# jobs chuvaMC na fila (recusa; FORCAR=1 compila assim mesmo). No fermi, no no de login (os nos de computacao nao tem
# cmake):
#     ./compilar.sh          (a primeira vez leva uns 15 a 30 min; depois, so o que mudou)
# Ambiente: BUILD (./build), JOBS (4), FORCAR=1.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
BUILD="${BUILD:-$AQUI/build}"
source "$NEOPZ/Projects/LearningMohrCoulomb2/cluster/modulos.sh"
if [ -z "${FORCAR:-}" ] && command -v squeue > /dev/null; then
    fila="$(squeue -h -u "${USER:-$(id -un)}" -n chuvaMC -o %i)" || { echo "o squeue nao respondeu: tente de novo (ou FORCAR=1)" >&2; exit 1; }
    [ -z "$fila" ] || { echo "ha jobs chuvaMC na fila usando a libpz.so de $BUILD: espere terminarem (ou FORCAR=1)" >&2; exit 1; }
fi
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON
      -DUSING_MKL=ON -DMKL_THREAD_MODEL=seq)
CATCH2="$NEOPZ/Projects/LearningMohrCoulomb2/cluster/deps/Catch2"
[ -f "$CATCH2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$CATCH2")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target SlopeRainMonteCarlo -j "${JOBS:-4}"
mkdir -p "$AQUI/bin"
cp "$BUILD/Projects/SlopeRandomFieldRefinement/SlopeRainMonteCarlo" "$AQUI/bin/.SlopeRainMonteCarlo.novo"
mv -f "$AQUI/bin/.SlopeRainMonteCarlo.novo" "$AQUI/bin/SlopeRainMonteCarlo" # troca atomica
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2> /dev/null || echo sem_git)"
[ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=no 2> /dev/null)" ] && v="$v+mod"
echo "${v}_$(date '+%F_%T')" > "$AQUI/bin/versao.txt"
echo "Compilado: $AQUI/bin/SlopeRainMonteCarlo ($(cat "$AQUI/bin/versao.txt"))"
