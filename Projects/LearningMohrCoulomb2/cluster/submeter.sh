#!/bin/bash
# submeter.sh [array]: submete a compilacao e, depois dela, o array de casos (padrao 1-21).
# Retorna em segundos; pode desligar a VPN logo em seguida. Os jobs seguem no cluster.
# Ambiente opcional: RECOMPILAR=1 (compila mesmo com bin/rahardjo pronto); EMAIL=<endereco> (o SLURM avisa por e-mail
# quando o array termina ou falha, se o cluster enviar e-mail). Os numeros dos jobs ficam em logs/jobs.txt (status.sh).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./modulos.sh
P=(); [ -n "${PARTICAO:-}" ] && P=(--partition="$PARTICAO")
MB=(); MA=()
if [ -n "${EMAIL:-}" ]; then MB=(--mail-user="$EMAIL" --mail-type=FAIL); MA=(--mail-user="$EMAIL" --mail-type=END,FAIL); fi
ARRAY="${1:-1-21}"
B=""
if [ -x bin/rahardjo ] && [ -z "${RECOMPILAR:-}" ]; then
  DEP=()
  echo "bin/rahardjo ja existe; pulando a compilacao (RECOMPILAR=1 para refazer)"
else
  # o CMake baixa o Catch2 ao configurar e os nos de computacao podem nao ter internet: baixa aqui, no no de login
  if [ ! -f deps/Catch2/CMakeLists.txt ]; then
    rm -rf deps/Catch2
    git -c advice.detachedHead=false clone -q --depth 1 --branch v2.13.10 https://github.com/catchorg/Catch2.git deps/Catch2
    echo "Catch2 v2.13.10 baixado em deps/Catch2"
  fi
  B=$(sbatch --parsable ${P[@]+"${P[@]}"} ${MB[@]+"${MB[@]}"} job_compilar.sbatch); B="${B%%;*}"
  DEP=(--dependency=afterok:"$B" --kill-on-invalid-dep=yes)
  echo "compilacao: job $B"
fi
A=$(sbatch --parsable ${P[@]+"${P[@]}"} ${MA[@]+"${MA[@]}"} ${DEP[@]+"${DEP[@]}"} --array="$ARRAY" job_slurm.sbatch); A="${A%%;*}"
echo "$(date '+%F %T') compilacao=${B:--} casos=$A array=$ARRAY" >> logs/jobs.txt
echo "casos: job array $A ($ARRAY)"
echo "acompanhe depois com: ./status.sh   (ou ./status.sh --slurm para ver tambem duracao e memoria no sacct)"
