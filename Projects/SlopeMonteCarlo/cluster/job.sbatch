#!/bin/bash
# job.sbatch: um no do Monte Carlo, PARALELO processos de 1 thread (um por nucleo fisico). O no NO de NOS roda os
# processos k = NO * PARALELO + p, p = 0 .. PARALELO - 1, de K = NOS * PARALELO: cada processo faz as tarefas i com
# i % K == k da fila (os deterministicos e as realizacoes PRIMEIRA..ULTIMA dos CASOS). No fim, junta os mc_resumo.csv.
# Submetido por submeter.sh, que fixa NO, NOS, PARALELO, PRIMEIRA, ULTIMA, CASOS e SAIDA.
#SBATCH --job-name=slopeMC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
# variaveis do ambiente do fermi, como no modelo de job do manual
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
: "${SAIDA:?}" "${NO:?}" "${NOS:?}" "${PARALELO:?}" "${PRIMEIRA:?}" "${ULTIMA:?}"
EXE="$PWD/bin/SlopeMonteCarloCluster"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$(hostname -s): $EXE ausente ou com bibliotecas faltando" >&2
    exit 1
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 SMC_VERSAO="$(cat bin/versao.txt 2>/dev/null || echo -)"
K=$((NOS * PARALELO))
mkdir -p "$SAIDA/logs"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): no $NO de $NOS, processos $((NO * PARALELO))..$((NO * PARALELO + PARALELO - 1)) de $K," \
    "realizacoes $PRIMEIRA..$ULTIMA, casos ${CASOS:-todos}, SAIDA=$SAIDA, versao $SMC_VERSAO"
for ((p = 0; p < PARALELO; p++)); do
    k=$((NO * PARALELO + p))
    "$EXE" rodar "$SAIDA" "$PRIMEIRA" "$ULTIMA" "$k" "$K" "${CASOS:-todos}" >> "$SAIDA/logs/p${k}de${K}.txt" 2>&1 &
done
wait
"$EXE" junta "$SAIDA"
echo "$(date '+%F %T') job $SLURM_JOB_ID: fim"
