#!/bin/bash
# compilar_hu2024.sh: compila so o SlopeRainMonteCarlo numa build e num bin proprios (build_hu/ e bin_hu/), para o
# hu2024_casos.sbatch nao depender da build/ e do bin/ das outras campanhas (que podem ter jobs rodando). Recusa compilar
# com jobs huA na fila (a libpz.so fica em build_hu/; FORCAR=1 compila assim mesmo). No no de login do fermi:
#     ./compilar_hu2024.sh      (a primeira vez leva uns 15 a 30 min)
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_hu}" BIN="$AQUI/bin_hu" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo \
    NOMES=huA exec "$AQUI/compilar.sh"
