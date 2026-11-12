#!/bin/bash
# job.sbatch: um no do Monte Carlo, PARALELO processos de 1 thread (um por nucleo fisico). O no NO de NOS roda os
# processos k = NO * PARALELO + p, p = 0 .. PARALELO - 1, de K = NOS * PARALELO: cada processo faz as tarefas i com
# i % K == k da fila (os deterministicos e as realizacoes PRIMEIRA..ULTIMA dos CASOS). No fim, junta os mc_resumo.csv.
# Submetido por submeter.sh, que fixa NO, NOS, PARALELO, PRIMEIRA, ULTIMA, CASOS, NREF e SAIDA.
# 15 min antes do limite de tempo, os processos terminam a tarefa corrente e param, e o job submete a sua continuacao
# (CONTINUAR=0 desliga): o mesmo NO de NOS, que comeca quando este acabar e faz o que falta. Se o SLURM a recusar (cota
# de jobs), fica o aviso $SAIDA/CONTINUACAO_PENDENTE, que status.sh mostra.
#SBATCH --job-name=slopeMC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=4-00:00:00
#SBATCH --signal=B:USR1@900
#SBATCH --output=logs/%x_%j.out
# o aviso do fim pode chegar antes de tudo estar pronto (um TEMPO curto): guardado aqui, tratado depois
trap 'AVISO_FIM=1' USR1
# variaveis do ambiente do fermi, como no modelo de job do manual
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
: "${SAIDA:?}" "${NO:?}" "${NOS:?}" "${PARALELO:?}" "${PRIMEIRA:?}" "${ULTIMA:?}" "${NREF:?}"
EXE="$PWD/bin/SlopeMonteCarloCluster"
if ! [ -x "$EXE" ] || ldd "$EXE" | grep -q "not found"; then
    echo "$(hostname -s): $EXE ausente ou com bibliotecas faltando" >&2
    exit 1
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 SMC_NREF="$NREF" SMC_VERSAO="$(cat bin/versao.txt 2>/dev/null || echo -)"
K=$((NOS * PARALELO))
mkdir -p "$SAIDA/logs"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): no $NO de $NOS, processos $((NO * PARALELO))..$((NO * PARALELO + PARALELO - 1)) de $K," \
    "realizacoes $PRIMEIRA..$ULTIMA, casos ${CASOS:-todos}, NREF=$NREF, SAIDA=$SAIDA, versao $SMC_VERSAO"
pids=()
for ((p = 0; p < PARALELO; p++)); do
    k=$((NO * PARALELO + p))
    "$EXE" rodar "$SAIDA" "$PRIMEIRA" "$ULTIMA" "$k" "$K" "${CASOS:-todos}" >> "$SAIDA/logs/p${k}de${K}.txt" 2>&1 &
    pids+=($!)
done

continuar() {
    trap '' USR1
    echo "$(date '+%F %T') limite de tempo em 15 min: os processos terminam a tarefa corrente e param"
    kill -USR1 "${pids[@]}" 2>/dev/null
    [ "${CONTINUAR:-1}" = 1 ] || return 0
    # o job avisado antes de 15 min de conta (um TEMPO curto demais) nao continua: a continuacao seria avisada do mesmo
    # jeito, e os jobs se repetiriam sem fazer nada
    if [ "$SECONDS" -lt 900 ]; then
        echo "$(date '+%F %T') sem continuacao: o aviso chegou com menos de 15 min de job (TEMPO curto demais)" >&2
        return 0
    fi
    local j tentativa
    # ate 3 tentativas, 2 min entre elas (a cota de jobs submetidos conta este job enquanto ele roda)
    for tentativa in 1 2 3; do
        if j="$(sbatch --parsable --dependency=afterany:"$SLURM_JOB_ID" ${TEMPO:+--time="$TEMPO"} \
            ${PARTICAO:+--partition="$PARTICAO"} job.sbatch 2> "logs/.sbatch_$SLURM_JOB_ID.err")"; then
            j="${j##*$'\n'}"
            echo "$(date '+%F %T') continuacao: job ${j%%;*}"
            echo "$(date '+%F %T') ${j%%;*} no $NO de $NOS: continuacao de $SLURM_JOB_ID" >> logs/jobs.txt
            rm -f "logs/.sbatch_$SLURM_JOB_ID.err"
            return 0
        fi
        j="$(cat "logs/.sbatch_$SLURM_JOB_ID.err" 2>/dev/null) $j"
        echo "$(date '+%F %T') continuacao recusada (tentativa $tentativa de 3): $j" >&2
        [ $tentativa -lt 3 ] && sleep 120
    done
    rm -f "logs/.sbatch_$SLURM_JOB_ID.err"
    echo "$(date '+%F %T') job $SLURM_JOB_ID (no $NO de $NOS): continuacao recusada pelo SLURM; rode ./submeter.sh $NOS" \
        "quando os jobs acabarem. Ultima resposta: $j" >> "$SAIDA/CONTINUACAO_PENDENTE"
}
vivos() {
    local p
    for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && return 0; done
    return 1
}
trap continuar USR1
[ -n "${AVISO_FIM:-}" ] && continuar
# wait volta quando chega o sinal (o trap roda); espera de novo ate os processos terminarem
while vivos; do wait; done
"$EXE" junta "$SAIDA"
echo "$(date '+%F %T') job $SLURM_JOB_ID: fim"
