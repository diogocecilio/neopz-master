#!/bin/bash
# submeter.sh [array]: submete a compilacao e, depois dela, o array de casos (padrao 1-21).
# Retorna em segundos; pode desligar a VPN logo em seguida. Os jobs seguem no cluster.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./modulos.sh
P=(); [ -n "${PARTICAO:-}" ] && P=(--partition="$PARTICAO")
ARRAY="${1:-1-21}"
if [ -x bin/rahardjo ] && [ -z "${RECOMPILAR:-}" ]; then
  DEP=()
  echo "bin/rahardjo ja existe; pulando a compilacao (RECOMPILAR=1 para refazer)"
else
  B=$(sbatch --parsable "${P[@]}" job_compilar.sbatch)
  DEP=(--dependency=afterok:"$B" --kill-on-invalid-dep=yes)
  echo "compilacao: job $B"
fi
A=$(sbatch --parsable "${P[@]}" "${DEP[@]}" --array="$ARRAY" job_slurm.sbatch)
echo "casos: job array $A ($ARRAY)"
echo "acompanhe depois com: squeue -u \$USER ; ls resultados/*/*/fs.csv"
