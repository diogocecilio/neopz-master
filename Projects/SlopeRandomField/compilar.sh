#!/bin/bash
# compilar.sh: compila a biblioteca NeoPZ e o alvo SlopeRandomField em ./build e copia o executavel para ./bin.
# Ambiente opcional: USAR_MKL=ON|OFF (padrao ON: PARDISO com ordenacao METIS), MKL_THREAD_MODEL (padrao seq), JOBS
# (padrao 4), BUILD (padrao ./build). Com os modulos/variaveis do MKL carregados (no fermi: source
# ../LearningMohrCoulomb2/cluster/modulos.sh, e compilar no no de login, que tem cmake).
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$PKG/../.." && pwd)"
BUILD="${BUILD:-$PKG/build}"
USAR_MKL="${USAR_MKL:-ON}"
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL="$USAR_MKL")
if [ "$USAR_MKL" = ON ]; then ARGS+=(-DMKL_THREAD_MODEL="${MKL_THREAD_MODEL:-seq}"); else ARGS+=(-DUSING_LAPACK=OFF); fi
# Catch2 local, se o pacote cluster de LearningMohrCoulomb2 ja o baixou (configurar sem internet)
CATCH2="$PKG/../LearningMohrCoulomb2/cluster/deps/Catch2"
[ -f "$CATCH2/CMakeLists.txt" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$CATCH2")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target SlopeRandomField -j "${JOBS:-4}"
# troca atomica: realizacoes rodando continuam com o executavel antigo
mkdir -p "$PKG/bin"
cp "$BUILD/Projects/SlopeRandomField/SlopeRandomField" "$PKG/bin/.SlopeRandomField.novo"
mv -f "$PKG/bin/.SlopeRandomField.novo" "$PKG/bin/SlopeRandomField"
echo "Compilado: $PKG/bin/SlopeRandomField"
