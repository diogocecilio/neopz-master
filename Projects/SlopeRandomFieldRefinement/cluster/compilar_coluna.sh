#!/bin/bash
# compilar_coluna.sh: compila o SlopeRainColumn (a coluna 1D de Cho 2014, coluna.cpp) e o SlopeRainMonteCarlo numa
# build e num bin proprios (build_coluna/ e bin_coluna/), para o cho2014_casos.sbatch (jobs cho14) nao depender das
# builds das outras campanhas (que podem ter jobs rodando). Recusa compilar com jobs cho14 ou liuCn na fila (a libpz.so
# fica em build_coluna/; FORCAR=1 compila assim mesmo). ALVOS=SlopeRandomFieldRefinement/SlopeRainColumn: so a coluna.
# No no de login do fermi:
#     ./compilar_coluna.sh      (a primeira vez leva uns 15 a 30 min)
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_coluna}" BIN="$AQUI/bin_coluna" \
    ALVOS="${ALVOS:-SlopeRandomFieldRefinement/SlopeRainColumn SlopeRandomFieldRefinement/SlopeRainMonteCarlo}" \
    NOMES=cho14,liuCn exec "$AQUI/compilar.sh"
