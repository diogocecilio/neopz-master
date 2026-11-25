#!/bin/bash
# compilar_det_chuva.sh: compila so o SlopeRainMonteCarlo para o det_chuva.sbatch (jobs detChuva), na build/ das
# campanhas (so o que mudou) e num bin proprio (bin_det/). Recusa com jobs que usam a build/ na fila (FORCAR=1 compila
# assim mesmo). No no de login do fermi:
#     ./compilar_det_chuva.sh
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$AQUI/bin_det" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo NOMES=detChuva,chuva2,chuvaMC,lsChuva2 \
    exec "$AQUI/compilar.sh"
