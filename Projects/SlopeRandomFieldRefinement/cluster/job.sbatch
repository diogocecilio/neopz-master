#!/bin/bash
# job.sbatch: um no do Monte Carlo: PARALELO processos de 1 thread (tarefas.sh), que pegam as tarefas livres da fila em
# $SAIDA. Submetido por submeter.sh, que fixa CAMPANHA (o executavel e os casos, campanha.sh), o nome do job, SAIDA,
# CASOS, IREF, BLOCO, ULTIMA e PARALELO; no primeiro job, FORM_M: antes da fila, os FORM que faltam (form.sh).
#SBATCH --job-name=chuvaMC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
# variaveis do ambiente do fermi, como no modelo de job do manual
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
source ./campanha.sh
: "${SAIDA:?}" "${PARALELO:?}" "${CASOS:?}" "${IREF:?}" "${BLOCO:?}"
export EXE="$PWD/bin/$PROG" CASOS_ARQ="$PWD/$CASOS_TXT"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$(hostname -s): $EXE ausente ou com bibliotecas faltando" >&2
    exit 1
fi
# um processo por nucleo, cada um com 1 thread (montagem e MKL)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA/logs"
echo "$(date '+%F %T') job $SLURM_JOB_ID ($CAMPANHA) em $(hostname -s): $PARALELO processos, casos $CASOS, IREF $IREF," \
    "blocos de $BLOCO${ULTIMA:+, ate a realizacao $ULTIMA}, SAIDA=$SAIDA, versao $(cat bin/versao.txt 2> /dev/null || echo -)"
[ -n "${FORM_M:-}" ] && ./form.sh
for ((p = 0; p < PARALELO; p++)); do
    ./tarefas.sh "${SLURM_JOB_ID}_$p" &
done
wait
echo "$(date '+%F %T') job $SLURM_JOB_ID: fim"
