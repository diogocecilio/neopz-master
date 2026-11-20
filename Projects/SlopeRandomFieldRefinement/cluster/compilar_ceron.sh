#!/bin/bash
# compilar_ceron.sh: compila so o SlopeDrawdownMonteCarlo numa build e num bin proprios (build_ceron/ e bin_ceron/), para
# o ls_ceron.sbatch nao depender da build/ e do bin/ das campanhas (que podem ter jobs rodando). Recusa compilar com
# jobs lsCeron na fila (a libpz.so fica em build_ceron/; FORCAR=1 compila assim mesmo). No no de login do fermi:
#     ./compilar_ceron.sh      (a primeira vez leva uns 15 a 30 min)
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_ceron}" BIN="$AQUI/bin_ceron" ALVOS=SlopeRandomFieldRefinement/SlopeDrawdownMonteCarlo \
    NOMES=lsCeron exec "$AQUI/compilar.sh"
