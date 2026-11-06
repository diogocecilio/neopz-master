#!/bin/bash
# job_mc.sbatch: um no do Monte Carlo de um projeto (PROJETO=V2_TaludesCho ou V4_TaludeLiu), um processo por nucleo
# fisico (com hyper-threading, dois por nucleo deixam cada um ~2x mais lento), rodando a fila de rodar_mc.sh: as
# realizacoes PRIMEIRA..ULTIMA de cada entrada de CASOS e, com DET=1, os deterministicos. Varios jobs iguais dividem a
# mesma fila sem repetir realizacoes. Submetido por submeter_mc.sh; a mao, de dentro de Projects/Cluster (logs/ e a
# pasta do pacote vem de la):
#   PROJETO=V2 CASOS="cho1 cho2" PRIMEIRA=1 ULTIMA=1000 DET=1 sbatch --job-name=verif-mc-V2 mc/job_mc.sbatch
# 15 min antes do limite de tempo, o job deixa de comecar realizacoes e, se ainda faltar alguma, submete a sua
# continuacao (CONTINUAR=0 desliga), que comeca quando ele acabar e retoma as interrompidas. Se o SLURM recusar a
# continuacao (cota de jobs), fica o aviso $SAIDA/<projeto>/CONTINUACAO_PENDENTE, que status_mc.sh mostra.
# No fermi: compilar antes no no de login (submeter_mc.sh faz isso).
#SBATCH --job-name=verif-mc
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
# a pasta Projects/Cluster: VERIF_CLUSTER (submeter_mc.sh a passa) ou a de onde o job foi submetido (o SLURM roda uma
# copia deste script em outra pasta)
PKG=""
for p in "${VERIF_CLUSTER:-}" "${SLURM_SUBMIT_DIR:-}" "${SLURM_SUBMIT_DIR:-}/.."; do
    if [ -n "$p" ] && [ -f "$p/comum.sh" ] && [ -f "$p/mc/rodar_mc.sh" ]; then PKG="$(cd "$p" && pwd)"; break; fi
done
[ -n "$PKG" ] || { echo "nao achei Projects/Cluster: submeta de la (ou com mc/submeter_mc.sh)" >&2; exit 1; }
export VERIF_CLUSTER="$PKG"
cd "$PKG" || exit 1
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
source "$PKG/mc/projetos.sh"
mc_projeto "${PROJETO:-}" || exit 1
export PROJETO="$MC_PROJETO"
# o nome verif-mc-<V2|V4> (o de status_mc.sh e da continuacao), tambem num job submetido a mao com o nome padrao
NOME="verif-mc-$MC_ID"
[ "${SLURM_JOB_NAME:-}" = "$NOME" ] || scontrol update JobId="$SLURM_JOB_ID" JobName="$NOME" > /dev/null 2>&1 || true
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
export SAIDA="${SAIDA:-$PKG/resultados_mc}"
D="$SAIDA/$MC_PROJETO"
mkdir -p "$D" "$PKG/logs"
falta="$(aviso_ldd "$PKG/bin/$MC_EXE")"
if [ -n "$falta" ]; then
    echo "$(hostname -s): $falta" >&2
    echo "$(date '+%F %T') $(hostname -s) job $SLURM_JOB_ID: saiu no inicio: $falta" >> "$D/AVISOS"
    exit 1
fi
export PARALELO="${PARALELO:-$(paralelo_job)}"
A="${PRIMEIRA:?defina PRIMEIRA}" B="${ULTIMA:?defina ULTIMA}"
echo "$(date '+%F %T') job $SLURM_JOB_ID em $(hostname -s): $MC_PROJETO, $PARALELO realizacoes ao mesmo tempo;" \
    "CASOS=${CASOS:-}, realizacoes $A..$B, DET=${DET:-0}, SAIDA=$SAIDA"

continuar() {
    trap '' USR1
    touch "$D/.parar.$SLURM_JOB_ID"
    echo "$(date '+%F %T') limite de tempo em 15 min: nenhuma realizacao nova neste job"
    [ "${CONTINUAR:-1}" = 1 ] || return 0
    local k j tentativa
    k="$("$PKG/mc/rodar_mc.sh" --faltam "$A" "$B")"
    [ "${k:-0}" -gt 0 ] || return 0
    # ate 3 tentativas, ESPERA_SBATCH_S (2 min) entre elas (a cota de jobs submetidos conta este job enquanto ele roda)
    for tentativa in 1 2 3; do
        if j="$(sbatch --parsable --dependency=afterany:"$SLURM_JOB_ID" --job-name="$NOME" \
            --output="$PKG/logs/%x_%j.out" ${TEMPO:+--time="$TEMPO"} ${PARTICAO:+--partition="$PARTICAO"} \
            "$PKG/mc/job_mc.sbatch" 2> "$PKG/logs/.sbatch_$SLURM_JOB_ID.err")"; then
            j="${j##*$'\n'}"
            echo "$(date '+%F %T') continuacao: job ${j%%;*} ($k realizacoes faltando)"
            echo "$(date '+%F %T') ${j%%;*} $MC_PROJETO continuacao de $SLURM_JOB_ID" >> "$PKG/logs/jobs_mc.txt"
            rm -f "$PKG/logs/.sbatch_$SLURM_JOB_ID.err"
            return 0
        fi
        j="$(cat "$PKG/logs/.sbatch_$SLURM_JOB_ID.err") $j"
        echo "$(date '+%F %T') continuacao recusada (tentativa $tentativa de 3): $j" >&2
        [ $tentativa -lt 3 ] && sleep "${ESPERA_SBATCH_S:-120}"
    done
    rm -f "$PKG/logs/.sbatch_$SLURM_JOB_ID.err"
    echo "$(date '+%F %T') job $SLURM_JOB_ID: continuacao recusada pelo SLURM com $k realizacoes faltando; rode" \
        "mc/submeter_mc.sh $MC_ID quando houver vaga na cota. Ultima resposta: $j" >> "$D/CONTINUACAO_PENDENTE"
}
trap continuar USR1
"$PKG/mc/rodar_mc.sh" "$A" "$B" &
pid=$!
[ -n "${AVISO_FIM:-}" ] && continuar
# wait volta quando chega o sinal (o trap roda); espera de novo ate rodar_mc.sh terminar
while kill -0 "$pid" 2>/dev/null; do wait "$pid"; done
rm -f "$D/.parar.$SLURM_JOB_ID"
echo "$(date '+%F %T') job $SLURM_JOB_ID: fim"
