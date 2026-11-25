#!/bin/bash
# compilar_ls_chuva2.sh: compila so o SlopeRainMonteCarlo numa build e num bin proprios (build_chuva2/ e bin_chuva2/), para
# o ls_chuva2.sbatch (jobs lsChuva2) nao depender da build/ e do bin/ das campanhas de Monte Carlo (que podem ter jobs
# rodando). Recusa compilar com jobs lsChuva2 na fila (a libpz.so fica em build_chuva2/; FORCAR=1 compila assim mesmo).
# No no de login do fermi:
#     ./compilar_ls_chuva2.sh      (a primeira vez leva uns 15 a 30 min)
# Para o LS usar o mesmo executavel do MC da chuva2 (bin/, do compilar.sh): EXE=$PWD/bin/SlopeRainMonteCarlo no submeter
# (cuidado: o compilar.sh nao espera os jobs lsChuva2).
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build_chuva2}" BIN="$AQUI/bin_chuva2" ALVOS=SlopeRandomFieldRefinement/SlopeRainMonteCarlo \
    NOMES=lsChuva2 exec "$AQUI/compilar.sh"
