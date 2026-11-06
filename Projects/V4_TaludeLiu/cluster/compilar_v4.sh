#!/bin/bash
# compilar_v4.sh: compila a biblioteca NeoPZ e os alvos V4_TaludeLiu e V4_testes em ./build e copia o executavel para
# ./bin (com a versao do git em bin/versao.txt, que vai para o resumo.csv de cada realizacao). Rode no no de login do
# cluster, com os modulos carregados (source ../../Cluster/modulos.sh). Ambiente opcional: USAR_MKL=ON|OFF (padrao ON:
# PARDISO), JOBS (padrao 4), BUILD (padrao ./build), CATCH2_DIR (copia local do Catch2 v2, se o no nao tiver internet),
# FORCAR=1 (compila mesmo com realizacoes rodando).
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEOPZ="$(cd "$PKG/../../.." && pwd)"
BUILD="${BUILD:-$PKG/build}"
USAR_MKL="${USAR_MKL:-ON}"
ARGS=(-S "$NEOPZ" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBUILD_PLASTICITY_MATERIALS=ON -DBUILD_PROJECTS=ON -DUSING_MKL="$USAR_MKL")
if [ "$USAR_MKL" = ON ]; then ARGS+=(-DMKL_THREAD_MODEL=seq); else ARGS+=(-DUSING_LAPACK=OFF); fi
[ -n "${CATCH2_DIR:-}" ] && ARGS+=(-DFETCHCONTENT_SOURCE_DIR_CATCH2="$CATCH2_DIR")
[ -f "$BUILD/CMakeCache.txt" ] || ARGS+=(-G "Unix Makefiles")
# o executavel usa a libpz.so de $BUILD: nao recompilar com realizacoes rodando
if [ -z "${FORCAR:-}" ] && pgrep -f "$PKG/bin/V4_TaludeLiu" > /dev/null 2>&1; then
    echo "ha realizacoes rodando com $PKG/bin/V4_TaludeLiu nesta maquina: espere terminarem (ou FORCAR=1)" >&2
    exit 1
fi
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target V4_TaludeLiu V4_testes -j "${JOBS:-4}"
mkdir -p "$PKG/bin"
cp "$BUILD/Projects/V4_TaludeLiu/V4_TaludeLiu" "$PKG/bin/.V4_TaludeLiu.novo"
mv -f "$PKG/bin/.V4_TaludeLiu.novo" "$PKG/bin/V4_TaludeLiu" # troca atomica
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2>/dev/null || echo sem_git)"
[ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=no 2>/dev/null)" ] && v="$v+mod"
echo "$v $(date '+%F_%T')" > "$PKG/bin/versao.txt"
echo "Compilado: $PKG/bin/V4_TaludeLiu ($(cat "$PKG/bin/versao.txt")); testes: $BUILD/Projects/V4_TaludeLiu/V4_testes"
