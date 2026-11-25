#!/bin/bash
# compilar_frag.sh: compila so o SlopeRainMonteCarlo para a campanha frag numa build e num bin proprios (build_frag/ e
# bin_frag/), sem tocar na build/ (a libpz.so das outras campanhas, ex.: o detChuva rodando). Recusa com jobs frag na
# fila (FORCAR=1 compila assim mesmo). No no de login do fermi (a primeira vez leva uns 15 a 30 min):
#     ./compilar_frag.sh
# Os jobs usam esse bin com BIN=bin_frag no submeter.sh: CAMPANHA=frag BIN=bin_frag PARALELO=24 ./submeter.sh 4
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_frag}" BIN="$AQUI/bin_frag" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo NOMES=frag \
    exec "$AQUI/compilar.sh"
