#!/bin/bash
# compilar_v1.sh: configura o NeoPZ (biblioteca + projetos de verificacao) e compila os alvos do V1 em $BUILD.
# Uso (no no de login do fermi, que tem cmake): source ../../Cluster/modulos.sh && ./compilar_v1.sh
# Ambiente opcional: BUILD (padrao <neopz>/build_verif), USAR_MKL=ON|OFF (padrao ON: PARDISO), JOBS (padrao 4).
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$AQUI/../../.." && pwd)"
BUILD="${BUILD:-$NEOPZ/build_verif}"
USAR_MKL="${USAR_MKL:-ON}"
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON
      -DUSING_MKL="$USAR_MKL")
if [ "$USAR_MKL" = ON ]; then ARGS+=(-DMKL_THREAD_MODEL="${MKL_THREAD_MODEL:-seq}"); else ARGS+=(-DUSING_LAPACK=OFF); fi
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target V1_LiakopoulosTrifasico V1_testes -j "${JOBS:-4}"
echo "Compilado: $BUILD/Projects/V1_LiakopoulosTrifasico/V1_LiakopoulosTrifasico"
