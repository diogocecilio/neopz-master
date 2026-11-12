#!/bin/bash
# submeter.sh [nos = 4]: submete <nos> jobs job.sbatch (um no inteiro cada) que dividem a fila de tarefas: os
# deterministicos e as realizacoes PRIMEIRA..ULTIMA (padrao 0..9999) dos CASOS (padrao todos: as Secoes 6.1 a 6.5), com
# NREF niveis de refinamento adaptativo (padrao 4). Volta em segundos: pode desconectar. Cada job vai ate o limite de
# tempo e submete sozinho a sua continuacao, ate a fila acabar. Rodar de novo, depois que os jobs acabarem, faz so o que
# falta (as tarefas feitas sao puladas); com ULTIMA maior, acrescenta realizacoes (ULTIMA=19999 ./submeter.sh 4).
# Ambiente: NREF, SAIDA (padrao $DATA/ceron_mc_h$NREF, na particao de dados do fermi; sem $DATA, ./resultados_h$NREF),
# PRIMEIRA, ULTIMA, CASOS (por exemplo CASOS=fase1 ou CASOS=ref,hw=0.9), PARALELO (processos por no: 12, os nucleos
# fisicos do fermi), TEMPO (limite do job, padrao 4-00:00:00, o maximo da particao short), PARTICAO, CONTINUAR=0 (sem
# continuacao no limite de tempo) e FORCAR=1 (submete mesmo com jobs slopeMC na fila).
# Cada SAIDA e de um NREF so: a pasta que ja tem resultados de outro NREF e recusada.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../LearningMohrCoulomb2/cluster/modulos.sh
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export NREF="${NREF:-4}"
SAIDA="${SAIDA:-${DATA:+$DATA/ceron_mc_h$NREF}}"
export SAIDA="${SAIDA:-$PWD/resultados_h$NREF}" PRIMEIRA="${PRIMEIRA:-0}" ULTIMA="${ULTIMA:-9999}" CASOS="${CASOS:-todos}"
export PARALELO="${PARALELO:-12}" NOS
[ -x bin/SlopeMonteCarloCluster ] || { echo "compile antes, no no de login: ./compilar.sh" >&2; exit 1; }
# a lista de casos e a configuracao da pasta de saida (senao os jobs sairiam logo no inicio)
SMC_NREF="$NREF" bin/SlopeMonteCarloCluster confere "$SAIDA" "$CASOS" || exit 2
# a divisao das tarefas depende de NOS e PARALELO: jobs de outra submissao ainda rodando repetiriam tarefas
fila="$(squeue -h -u "${USER:-$(id -un)}" -n slopeMC -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -n "$fila" ] && [ -z "${FORCAR:-}" ]; then
    echo "ha $(echo "$fila" | wc -l) job(s) slopeMC na fila: espere terminarem (ou scancel -n slopeMC) e submeta de novo" >&2
    exit 1
fi
mkdir -p logs "$SAIDA"
rm -f "$SAIDA/PARAR"
ok=0
for ((no = 0; no < NOS; no++)); do
    j="$(NO=$no sbatch --parsable ${PARTICAO:+--partition="$PARTICAO"} ${TEMPO:+--time="$TEMPO"} job.sbatch)" || {
        echo "o SLURM recusou o job do no $no (cota?): os ja submetidos seguem; submeta o resto depois" >&2
        break
    }
    ok=$((ok + 1))
    echo "job ${j%%;*}: no $no de $NOS"
    echo "$(date '+%F %T') ${j%%;*} no $no de $NOS PARALELO=$PARALELO NREF=$NREF realizacoes $PRIMEIRA..$ULTIMA CASOS=$CASOS SAIDA=$SAIDA" >> logs/jobs.txt
done
[ "$ok" -gt 0 ] && rm -f "$SAIDA/CONTINUACAO_PENDENTE"
echo "resultados em $SAIDA; acompanhe com ./status.sh"
