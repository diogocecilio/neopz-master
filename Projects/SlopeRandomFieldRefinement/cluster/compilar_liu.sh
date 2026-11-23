#!/bin/bash
# compilar_liu.sh: compila so o SlopeRainMonteCarlo numa build e num bin proprios (build_liu/ e bin_liu/), para os casos
# B (liu2024_casos.sbatch, jobs liuB) e C (liuwang2023_casos.sbatch, jobs liuC), sem tocar na build_hu/ do caso A nem na
# build/ das outras campanhas (que podem ter jobs rodando). Recusa compilar com jobs liuB ou liuC na fila (a libpz.so
# fica em build_liu/; FORCAR=1 compila assim mesmo). No no de login do fermi:
#     ./compilar_liu.sh      (a primeira vez leva uns 15 a 30 min)
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_liu}" BIN="$AQUI/bin_liu" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo \
    NOMES=liuB,liuC exec "$AQUI/compilar.sh"
