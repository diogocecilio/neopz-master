#!/bin/bash
# submeter_casoC.sh: os 4 jobs do diagnostico do caso C (casoC_diag.sbatch), um por grupo de variantes:
#   tri3 (ref, kw4, kw3, kw2: 8 rodadas, 1 thread), tri4 (tri4, tri4kw4: 4 rodadas, 6 threads),
#   camadas (cam, camkw4: 4 rodadas, 6 threads), tri5 (tri5, tri5kw4: 4 rodadas, 6 threads).
# Resumo: ./casoC_diag.sbatch resumo
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
mkdir -p logs
VARS=ref,kw4,kw3,kw2 PZT=1 sbatch --export=ALL casoC_diag.sbatch
VARS=tri4,tri4kw4 PZT=6 sbatch --export=ALL casoC_diag.sbatch
VARS=cam,camkw4 PZT=6 sbatch --export=ALL casoC_diag.sbatch
VARS=tri5,tri5kw4 PZT=6 sbatch --export=ALL casoC_diag.sbatch
