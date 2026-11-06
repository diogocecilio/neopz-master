#!/bin/bash
# job_v4.sbatch: um no do Monte Carlo do V4 (SLURM), um processo por nucleo fisico, rodando rodar_v4.sh nas realizacoes
# PRIMEIRA..ULTIMA. Compile antes no no de login (./compilar_v4.sh) e crie a pasta logs/. Exemplos (na pasta cluster/):
#   mkdir -p logs; CASO=s1 ANGULO=0 PRIMEIRA=1 ULTIMA=1000 DET=1 sbatch job_v4.sbatch
#   CASO=s2 ANGULO=0 PRIMEIRA=1 ULTIMA=1000 sbatch job_v4.sbatch            # roda o det antes e reusa o permanente
#   CASO=s3 ANGULO=90 RHO=-0.5 PRIMEIRA=1 ULTIMA=1000 sbatch job_v4.sbatch  # com a correlacao cruzada do artigo
# Varios jobs: de faixas disjuntas (PRIMEIRA..ULTIMA), ou da mesma faixa um depois do outro (retoma as que faltam).
# Tempo por realizacao (malha de 0,5 m, 1 CPU): ~1,5 min (s1, s3: transiente + SRM) e ~1 min (s2, permanente lido).
#SBATCH --job-name=v4liu
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
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): CASO=${CASO:-s3} ANGULO=${ANGULO:-0} RHO=${RHO:-0}" \
    "realizacoes ${PRIMEIRA:?defina PRIMEIRA}..${ULTIMA:?defina ULTIMA}, $PARALELO ao mesmo tempo"
./rodar_v4.sh "$PRIMEIRA" "$ULTIMA"
