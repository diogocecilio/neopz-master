#!/bin/bash
# job_v2.sbatch: um no do Monte Carlo do V2 (SLURM), um processo por nucleo fisico, rodando rodar_v2.sh nas realizacoes
# PRIMEIRA..ULTIMA. Compile antes no no de login (./compilar_v2.sh) e crie a pasta logs/. Exemplos (na pasta cluster/):
#   mkdir -p logs; CASO=cho1 NUCLEO=gauss PRIMEIRA=1 ULTIMA=1000 DET=1 sbatch job_v2.sbatch
#   CASO=cho1 NUCLEO=exp PRIMEIRA=1 ULTIMA=1000 sbatch job_v2.sbatch
# Varios jobs: de faixas disjuntas (PRIMEIRA..ULTIMA), ou da mesma faixa um depois do outro (retoma as que faltam).
#SBATCH --job-name=v2cho
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../Cluster/modulos.sh
TPC="$(LC_ALL=C lscpu | awk -F: '/^Thread\(s\) per core/ { gsub(/ /, "", $2); print $2 }')"
case "$TPC" in '' | *[!0-9]* | 0) TPC=1 ;; esac
export PARALELO="${PARALELO:-$((SLURM_CPUS_PER_TASK / TPC))}"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): CASO=${CASO:-cho1} NUCLEO=${NUCLEO:-gauss}" \
    "realizacoes ${PRIMEIRA:?defina PRIMEIRA}..${ULTIMA:?defina ULTIMA}, $PARALELO ao mesmo tempo"
./rodar_v2.sh "$PRIMEIRA" "$ULTIMA"
