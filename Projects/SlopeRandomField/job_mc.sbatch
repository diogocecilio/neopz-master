#!/bin/bash
# job_mc.sbatch: realizacoes PRIMEIRA..ULTIMA de uma combinacao, num no, um processo por nucleo fisico (com hyper-
# threading, dois por nucleo deixam cada um ~2x mais lento). A combinacao e a estatistica vem do ambiente, por exemplo:
#   SOLO=f50 I_MM_H=36 PRIMEIRA=0 ULTIMA=99 sbatch job_mc.sbatch
# No fermi: compilar antes no no de login (./compilar.sh, com ../LearningMohrCoulomb2/cluster/modulos.sh carregado).
#SBATCH --job-name=slopeRF
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../LearningMohrCoulomb2/cluster/modulos.sh
if ! [ -x bin/SlopeRandomField ] || ldd bin/SlopeRandomField | grep -q "not found"; then
  echo "$(hostname -s): bin/SlopeRandomField ausente ou com bibliotecas faltando:" >&2
  [ -x bin/SlopeRandomField ] && ldd bin/SlopeRandomField | grep "not found" >&2
  exit 1
fi
TPC="$(lscpu | awk -F: '/Thread\(s\) per core/ { gsub(/ /, "", $2); print $2 }')"
export PARALELO="${PARALELO:-$((SLURM_CPUS_PER_TASK / ${TPC:-1}))}"
./rodar_mc.sh "${PRIMEIRA:?defina PRIMEIRA}" "${ULTIMA:?defina ULTIMA}"
