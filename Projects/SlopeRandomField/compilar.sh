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
# o executavel em bin/ carrega a libpz.so de $BUILD: recompilar com realizacoes rodando (aqui ou em jobs) as deixaria com
# a biblioteca nova e o executavel velho. FORCAR=1 compila mesmo assim
if [ -z "${FORCAR:-}" ]; then
    if pgrep -f "$PKG/bin/SlopeRandomField" > /dev/null 2>&1; then
        echo "ha realizacoes rodando nesta maquina com $PKG/bin/SlopeRandomField: espere terminarem (ou FORCAR=1)" >&2
        exit 1
    fi
    if command -v squeue > /dev/null; then
        if ! fila="$(squeue -h -u "${USER:-$(id -un)}" -n slopeRF -o %i 2>&1)"; then
            echo "o squeue nao respondeu ($fila): sem saber se ha jobs usando a libpz.so de $BUILD, nao compila (ou FORCAR=1)" >&2
            exit 1
        fi
        if [ -n "$fila" ]; then
            echo "ha jobs slopeRF no SLURM usando a libpz.so de $BUILD: espere terminarem (ou FORCAR=1)" >&2
            exit 1
        fi
    fi
fi
cmake "${ARGS[@]}"
cmake --build "$BUILD" --target SlopeRandomField -j "${JOBS:-4}"
# troca atomica: realizacoes rodando continuam com o executavel antigo
mkdir -p "$PKG/bin"
cp "$BUILD/Projects/SlopeRandomField/SlopeRandomField" "$PKG/bin/.SlopeRandomField.novo"
mv -f "$PKG/bin/.SlopeRandomField.novo" "$PKG/bin/SlopeRandomField"
# a versao: o commit do git (com "+mod" se havia alteracoes locais) e a data; vai para o resumo.csv de cada realizacao
v="$(git -C "$NEOPZ" rev-parse --short HEAD 2>/dev/null || echo sem_git)"
[ -n "$(git -C "$NEOPZ" status --porcelain --untracked-files=no 2>/dev/null)" ] && v="$v+mod"
echo "$v $(date '+%F_%T')" > "$PKG/bin/versao.txt"
echo "Compilado: $PKG/bin/SlopeRandomField"
