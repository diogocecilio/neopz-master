#!/bin/bash
# compilar.sh: compila a NeoPZ, o SlopeRainMonteCarlo (chuva) e o SlopeDrawdownMonteCarlo (rebaixamento do Ceron) numa
# pasta de build PROPRIA (padrao: ./build), sem tocar na build/ do SlopeMonteCarlo, e copia os executaveis para ./bin
# com a versao (o commit do git, "+mod" se havia alteracoes locais, e a data) em bin/versao.txt. Os jobs usam bin/ (a
# troca e atomica), mas a libpz.so continua em build/: nao compile com jobs chuvaMC ou ceronMC na fila (recusa;
# FORCAR=1 compila assim mesmo). No fermi, no no de login (os nos de computacao nao tem
# cmake):
#     ./compilar.sh          (a primeira vez leva uns 15 a 30 min; depois, so o que mudou)
# Ambiente: BUILD (./build), JOBS (4), FORCAR=1.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
BUILD="${BUILD:-$AQUI/build}"
source "$NEOPZ/Projects/LearningMohrCoulomb2/cluster/modulos.sh"
source "$AQUI/campanha.sh"
if [ -z "${FORCAR:-}" ] && command -v squeue > /dev/null; then
    fila="$(squeue -h -u "${USER:-$(id -un)}" -n "$NOMES_TODOS" -o %i)" || { echo "o squeue nao respondeu: tente de novo (ou FORCAR=1)" >&2; exit 1; }
    [ -z "$fila" ] || { echo "ha jobs ($NOMES_TODOS) na fila usando a libpz.so de $BUILD: espere terminarem (ou FORCAR=1)" >&2; exit 1; }
fi
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON
      -DUSING_MKL=ON -DMKL_THREAD_MODEL=seq)
CATCH2="$NEOPZ/Projects/LearningMohrCoulomb2/cluster/deps/Catch2"
[ -f "$CATCH2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$CATCH2")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
cmake "${ARGS[@]}"
mkdir -p "$AQUI/bin"
for prog in SlopeRainMonteCarlo SlopeDrawdownMonteCarlo Cho2010; do
    cmake --build "$BUILD" --target "$prog" -j "${JOBS:-4}"
    exe="$(find "$BUILD/Projects" -maxdepth 2 -type f -name "$prog" -perm -u+x | head -n 1)"
    [ -n "$exe" ] || { echo "nao achei o executavel $prog em $BUILD/Projects" >&2; exit 1; }
    cp "$exe" "$AQUI/bin/.$prog.novo"
    mv -f "$AQUI/bin/.$prog.novo" "$AQUI/bin/$prog" # troca atomica
done
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2> /dev/null || echo sem_git)"
[ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=no 2> /dev/null)" ] && v="$v+mod"
echo "${v}_$(date '+%F_%T')" > "$AQUI/bin/versao.txt"
echo "Compilado: $AQUI/bin/SlopeRainMonteCarlo, SlopeDrawdownMonteCarlo e Cho2010 ($(cat "$AQUI/bin/versao.txt"))"
