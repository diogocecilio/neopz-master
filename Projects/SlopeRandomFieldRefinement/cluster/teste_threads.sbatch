#!/bin/bash
# teste_threads.sbatch: no mesmo no, roda 12 realizacoes do Monte Carlo com chuva em paralelo (uma por nucleo), primeiro
# com a montagem serial (PZ_THREADS=0) e depois com 8 threads por processo (PZ_THREADS=8), e compara os tempos.
# Submeter da pasta cluster/:   sbatch teste_threads.sbatch
#SBATCH --job-name=testeThreads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --exclusive
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=teste_threads_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
EXE="$PWD/build/Projects/SlopeRandomFieldRefinement/SlopeRainMonteCarlo"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
N=${N:-12}
SAIDA="${SAIDA:-$PWD/teste_threads_$SLURM_JOB_ID}"
echo "$(date '+%F %T') no $(hostname -s), $N processos, saida $SAIDA"
for T in 0 8; do
    t0=$(date +%s)
    for ((p = 0; p < N; p++)); do
        d="$SAIDA/T$T/p$p"; mkdir -p "$d"
        ( cd "$d" && PZ_THREADS=$T /usr/bin/time -f "%e" -o tempo.txt "$EXE" n=1 primeira=$p transiente=tri2 uT=2 \
          fator=tri1 nref=4 > out.txt 2>&1 ) &
    done
    wait
    t1=$(date +%s)
    echo "PZ_THREADS=$T: $N realizacoes em $((t1 - t0)) s (parede); por realizacao:" \
         "$(cat "$SAIDA"/T$T/p*/tempo.txt | awk '{s+=$1} END {printf "%.1f s em media", s/NR}')"
    grep -h "realizacao  " "$SAIDA"/T$T/p*/out.txt | tr -d '*'
done
echo "$(date '+%F %T') fim"
