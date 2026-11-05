#!/bin/bash
# job_mc.sbatch: um no do Monte Carlo, um processo por nucleo fisico (com hyper-threading, dois por nucleo deixam cada um
# ~2x mais lento), rodando a fila de rodar_mc.sh: as realizacoes PRIMEIRA..ULTIMA de cada combinacao de COMBS e, com
# DET=1, os deterministicos. Varios jobs iguais dividem a mesma fila sem repetir realizacoes. Submetido por
# submeter_mc.sh; a mao: COMBS="f10:9 f50:80" PRIMEIRA=0 ULTIMA=999 sbatch job_mc.sbatch.
# 15 min antes do limite de tempo, o job deixa de comecar realizacoes e, se ainda faltar alguma, submete a sua
# continuacao (CONTINUAR=0 desliga), que comeca quando ele acabar e retoma as interrompidas dos seus checkpoints. Se o
# SLURM recusar a continuacao (cota de jobs), fica o aviso $SAIDA/CONTINUACAO_PENDENTE, que status_mc.sh mostra.
# No fermi: compilar antes no no de login (submeter_mc.sh faz isso).
#SBATCH --job-name=slopeRF
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=4-00:00:00
#SBATCH --signal=B:USR1@900
#SBATCH --output=logs/%x_%j.out
# o aviso do fim pode chegar antes de tudo estar pronto (um TEMPO curto): guardado aqui, tratado depois
trap 'AVISO_FIM=1' USR1
# variaveis adicionais do ambiente do fermi, como no modelo de job do manual (antes do cd: o script pode mudar de pasta)
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../LearningMohrCoulomb2/cluster/modulos.sh
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
export SAIDA="${SAIDA:-$PWD/resultados_mc}"
mkdir -p "$SAIDA"
if ! [ -x bin/SlopeRandomField ] || ldd bin/SlopeRandomField | grep -q "not found"; then
    falta="$( { [ -x bin/SlopeRandomField ] && ldd bin/SlopeRandomField | grep "not found" | tr '\n' ' '; } || echo "sem o executavel")"
    echo "$(hostname -s): bin/SlopeRandomField ausente ou com bibliotecas faltando: $falta" >&2
    echo "$(date '+%F %T') $(hostname -s) job $SLURM_JOB_ID: saiu no inicio, bin/SlopeRandomField ausente ou com bibliotecas faltando: $falta" >> "$SAIDA/AVISOS"
    exit 1
fi
TPC="$(LC_ALL=C lscpu | awk -F: '/^Thread\(s\) per core/ { gsub(/ /, "", $2); print $2 }')"
case "$TPC" in '' | *[!0-9]* | 0) TPC=1 ;; esac
export PARALELO="${PARALELO:-$((SLURM_CPUS_PER_TASK / TPC))}"
A="${PRIMEIRA:?defina PRIMEIRA}" B="${ULTIMA:?defina ULTIMA}"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): $PARALELO realizacoes ao mesmo tempo; COMBS=${COMBS:-}" \
    "realizacoes $A..$B, DET=${DET:-0}, SAIDA=$SAIDA"

continuar() {
    trap '' USR1
    touch "$SAIDA/.parar.$SLURM_JOB_ID"
    echo "$(date '+%F %T') limite de tempo em 15 min: nenhuma realizacao nova neste job"
    [ "${CONTINUAR:-1}" = 1 ] || return 0
    local k j tentativa
    k="$(./rodar_mc.sh --faltam "$A" "$B")"
    [ "${k:-0}" -gt 0 ] || return 0
    # ate 3 tentativas, 2 min entre elas (a cota de jobs submetidos conta este job enquanto ele roda)
    for tentativa in 1 2 3; do
        if j="$(sbatch --parsable --dependency=afterany:"$SLURM_JOB_ID" ${TEMPO:+--time="$TEMPO"} \
            ${PARTICAO:+--partition="$PARTICAO"} job_mc.sbatch 2> "logs/.sbatch_$SLURM_JOB_ID.err")"; then
            j="${j##*$'\n'}"
            echo "$(date '+%F %T') continuacao: job ${j%%;*} ($k realizacoes faltando)"
            echo "$(date '+%F %T') ${j%%;*} continuacao de $SLURM_JOB_ID" >> logs/jobs_mc.txt
            rm -f "logs/.sbatch_$SLURM_JOB_ID.err"
            return 0
        fi
        j="$(cat "logs/.sbatch_$SLURM_JOB_ID.err") $j"
        echo "$(date '+%F %T') continuacao recusada (tentativa $tentativa de 3): $j" >&2
        [ $tentativa -lt 3 ] && sleep 120
    done
    echo "$(date '+%F %T') job $SLURM_JOB_ID: continuacao recusada pelo SLURM com $k realizacoes faltando; rode" \
        "./submeter_mc.sh quando houver vaga na cota. Ultima resposta: $j" >> "$SAIDA/CONTINUACAO_PENDENTE"
}
trap continuar USR1
./rodar_mc.sh "$A" "$B" &
pid=$!
[ -n "${AVISO_FIM:-}" ] && continuar
# wait volta quando chega o sinal (o trap roda); espera de novo ate rodar_mc.sh terminar
while kill -0 "$pid" 2>/dev/null; do wait "$pid"; done
rm -f "$SAIDA/.parar.$SLURM_JOB_ID"
echo "$(date '+%F %T') job $SLURM_JOB_ID: fim"
