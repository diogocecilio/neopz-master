#!/bin/bash
# Compila a biblioteca NeoPZ deste repositorio e o driver rahardjo.cpp, com as flags e a ligacao do projeto
# LearningMohrCoulomb2. Rode no no de login, com os modulos carregados. Ambiente opcional:
#   USAR_MKL=ON|OFF (padrao ON; OFF usa a skyline, mais lenta)   MKL_ROOT=<pasta do MKL>
#   JOBS=<nucleos, padrao 4>   BUILD=<pasta de build, padrao ./build>
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$PKG/../../.." && pwd)"
BUILD="${BUILD:-$PKG/build}"
USAR_MKL="${USAR_MKL:-ON}"
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL="$USAR_MKL")
[ -n "${MKL_ROOT:-}" ] && ARGS+=(-DMKL_ROOT="$MKL_ROOT")
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
