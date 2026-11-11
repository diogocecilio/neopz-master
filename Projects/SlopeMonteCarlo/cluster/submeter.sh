#!/bin/bash
# submeter.sh [nos = 4]: submete <nos> jobs job.sbatch (um no inteiro cada) que dividem a fila de tarefas: os
# deterministicos e as realizacoes PRIMEIRA..ULTIMA (padrao 0..999) dos CASOS (padrao todos: as Secoes 6.1 a 6.5). Volta
# em segundos: pode desconectar. Rodar de novo, depois que os jobs acabarem, faz so o que falta (as tarefas feitas sao
# puladas); com ULTIMA maior, acrescenta realizacoes (ULTIMA=1999 ./submeter.sh 4).
# Ambiente: SAIDA (padrao $DATA/ceron_mc, a particao de dados do fermi; sem $DATA, ./resultados), PRIMEIRA, ULTIMA,
# CASOS (por exemplo CASOS=fase1 ou CASOS=ref,hw=0.9), PARALELO (processos por no: 12, os nucleos fisicos do fermi),
# TEMPO (limite do job, padrao 2-00:00:00), PARTICAO e FORCAR=1 (submete mesmo com jobs slopeMC na fila).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../LearningMohrCoulomb2/cluster/modulos.sh
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
SAIDA="${SAIDA:-${DATA:+$DATA/ceron_mc}}"
export SAIDA="${SAIDA:-$PWD/resultados}" PRIMEIRA="${PRIMEIRA:-0}" ULTIMA="${ULTIMA:-999}" CASOS="${CASOS:-todos}"
export PARALELO="${PARALELO:-12}" NOS
[ -x bin/SlopeMonteCarloCluster ] || { echo "compile antes, no no de login: ./compilar.sh" >&2; exit 1; }
bin/SlopeMonteCarloCluster casos "$CASOS" > /dev/null || exit 2 # confere a lista de casos
# a divisao das tarefas depende de NOS e PARALELO: jobs de outra submissao ainda rodando repetiriam tarefas
fila="$(squeue -h -u "${USER:-$(id -un)}" -n slopeMC -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -n "$fila" ] && [ -z "${FORCAR:-}" ]; then
    echo "ha $(echo "$fila" | wc -l) job(s) slopeMC na fila: espere terminarem (ou scancel -n slopeMC) e submeta de novo" >&2
    exit 1
fi
mkdir -p logs "$SAIDA"
rm -f "$SAIDA/PARAR"
for ((no = 0; no < NOS; no++)); do
    j="$(NO=$no sbatch --parsable ${PARTICAO:+--partition="$PARTICAO"} ${TEMPO:+--time="$TEMPO"} job.sbatch)" || {
        echo "o SLURM recusou o job do no $no (cota?): os ja submetidos seguem; submeta o resto depois" >&2
        break
    }
    echo "job ${j%%;*}: no $no de $NOS"
    echo "$(date '+%F %T') ${j%%;*} no $no de $NOS PARALELO=$PARALELO realizacoes $PRIMEIRA..$ULTIMA CASOS=$CASOS SAIDA=$SAIDA" >> logs/jobs.txt
done
echo "resultados em $SAIDA; acompanhe com ./status.sh"
