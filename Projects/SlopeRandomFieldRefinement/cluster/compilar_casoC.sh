#!/bin/bash
# compilar_casoC.sh: compila so o SlopeRainMonteCarlo para o casoC_diag.sbatch (jobs casoCd), na build_hu/ do caso A
# (a campanha huA ja terminou: so recompila o que mudou) e num bin proprio (bin_casoC/), sem tocar na build_liu/ dos
# jobs liuB que estao rodando. Recusa com jobs huA ou casoCd na fila (FORCAR=1 compila assim mesmo). No no de login:
#     ./compilar_casoC.sh
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_hu}" BIN="$AQUI/bin_casoC" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo \
    NOMES=huA,casoCd exec "$AQUI/compilar.sh"
