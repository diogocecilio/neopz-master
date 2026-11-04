#!/bin/bash
# Compila a biblioteca NeoPZ deste repositorio e o driver rahardjo.cpp, com as flags e a ligacao do projeto
# LearningMohrCoulomb2. Roda no job de compilacao (job_compilar.sbatch) ou direto, com os modulos carregados. Ambiente opcional:
#   USAR_MKL=ON|OFF (padrao ON; OFF usa a skyline, mais lenta)   MKL_ROOT=<pasta do MKL>
#   JOBS=<nucleos, padrao 4>   BUILD=<pasta de build, padrao ./build>
#   MKL_THREAD_MODEL=seq|omp|tbb (padrao seq). Cada caso roda num processo de 1 CPU, e o MKL sequencial nao abre threads
#   nem depende do TBB. Sem isso, o NeoPZ escolhe tbb sempre que acha o TBB (no fermi o modulo tbb e pre-requisito do mkl)
#   e compila a biblioteca inteira com USING_TBB.
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$PKG/../../.." && pwd)"
BUILD="${BUILD:-$PKG/build}"
USAR_MKL="${USAR_MKL:-ON}"
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL="$USAR_MKL")
[ -n "${MKL_ROOT:-}" ] && ARGS+=(-DMKL_ROOT="$MKL_ROOT")
[ "$USAR_MKL" = ON ] && ARGS+=(-DMKL_THREAD_MODEL="${MKL_THREAD_MODEL:-seq}")
# Projects/CMakeLists.txt baixa o Catch2 do GitHub ao configurar, e os nos de computacao podem nao ter internet:
# submeter.sh o baixa antes, no no de login, em deps/Catch2, e o CMake usa essa copia sem acessar a rede
[ -f "$PKG/deps/Catch2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$PKG/deps/Catch2")
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target LearningMohrCoulomb2 -j "${JOBS:-4}"
# o driver usa as flags e a linha de ligacao geradas para o projeto
P="$BUILD/Projects/LearningMohrCoulomb2"
D="$P/CMakeFiles/LearningMohrCoulomb2.dir"
DEFS="$(grep '^CXX_DEFINES' "$D/flags.make" | sed 's/^CXX_DEFINES = //')"
INCS="$(grep '^CXX_INCLUDES' "$D/flags.make" | sed 's/^CXX_INCLUDES = //')"
LIBS="$(sed -E 's/.* -o LearningMohrCoulomb2 //' "$D/link.txt")"
mkdir -p "$PKG/bin"
cd "$P"
eval "${CXX:-c++}" $DEFS $INCS -O3 -DNDEBUG -std=gnu++17 '"$PKG/fonte/rahardjo.cpp"' -o '"$PKG/bin/rahardjo"' $LIBS
echo "Compilado: $PKG/bin/rahardjo"
