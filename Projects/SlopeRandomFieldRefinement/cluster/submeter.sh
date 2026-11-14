#!/bin/bash
# submeter.sh [nos = 4]: submete <nos> jobs job.sbatch (um no inteiro cada), que dividem a fila de tarefas dos CASOS de
# casos.txt. Volta em segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam as tarefas
# livres. Sem jobs chuvaMC na fila, libera antes as tarefas interrompidas (as travas sem mapa), que sao refeitas.
# Ambiente:
#   SAIDA     pasta dos resultados (padrao $DATA/chuva_mc, a particao de dados do fermi; sem $DATA, ./resultados)
#   CASOS     todos (padrao), ids ou fases separados por virgula: CASOS=fase1, CASOS=ref,I=36
#   IREF      intensidade de referencia em mm/h (padrao 31.05); fica gravada em $SAIDA/campanha.txt e nao pode mudar
#   BLOCO     realizacoes por tarefa (padrao 20, umas 2 h no fermi); tambem fixo por SAIDA
#   ULTIMA    ultima realizacao de todos os casos escolhidos (ex.: ULTIMA=9999 CASOS=ref acrescenta realizacoes)
#   PARALELO  processos por no (padrao 12, os nucleos fisicos; 24 usa o hyperthreading)
#   TEMPO     limite de cada job (padrao 2-00:00:00), PARTICAO
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../LearningMohrCoulomb2/cluster/modulos.sh
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
SAIDA="${SAIDA:-${DATA:+$DATA/chuva_mc}}"
export SAIDA="${SAIDA:-$PWD/resultados}" CASOS="${CASOS:-todos}" IREF="${IREF:-31.05}" BLOCO="${BLOCO:-20}"
export PARALELO="${PARALELO:-12}"
[ -n "${ULTIMA:-}" ] && export ULTIMA
[ -x bin/SlopeRainMonteCarlo ] || { echo "compile antes, no no de login: ./compilar.sh" >&2; exit 1; }

# os casos pedidos existem?
if [ "$CASOS" != todos ]; then
    for c in ${CASOS//,/ }; do
        [[ "$c" == fase[123] ]] && continue
        awk -v c="$c" '$1 == c { achou = 1 } END { exit !achou }' casos.txt || { echo "caso desconhecido: $c" >&2; exit 2; }
    done
fi

# IREF e BLOCO definem os arquivos de cada caso: fixos durante a campanha
mkdir -p logs "$SAIDA"
if [ -f "$SAIDA/campanha.txt" ]; then
    read -r i0 b0 < "$SAIDA/campanha.txt"
    if [ "$i0" != "$IREF" ] || [ "$b0" != "$BLOCO" ]; then
        echo "$SAIDA e da campanha IREF=$i0 BLOCO=$b0: use os mesmos valores (ou outra SAIDA)" >&2
        exit 1
    fi
else
    echo "$IREF $BLOCO" > "$SAIDA/campanha.txt"
fi

fila="$(squeue -h -u "${USER:-$(id -un)}" -n chuvaMC -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    # nenhum job na fila: uma trava sem o mapa e de uma tarefa interrompida (fim do tempo, queda do no): libera
    n=0
    while IFS= read -r -d '' trava; do
        t="$(basename "$trava")" caso="$(dirname "$(dirname "$trava")")"
        [ -f "$caso/mc_mapa_$t.csv" ] || { rm -rf "$trava"; n=$((n + 1)); }
    done < <(find "$SAIDA" -path '*/.tarefas/*' -mindepth 1 -print0 2> /dev/null)
    [ "$n" -eq 0 ] || echo "liberadas $n tarefas interrompidas (serao refeitas)"
else
    echo "ha $(echo "$fila" | wc -l) job(s) chuvaMC na fila: os novos entram na mesma fila de tarefas"
fi
rm -f "$SAIDA/PARAR"
for ((no = 0; no < NOS; no++)); do
    j="$(sbatch --parsable ${PARTICAO:+--partition="$PARTICAO"} ${TEMPO:+--time="$TEMPO"} job.sbatch)" || {
        echo "o SLURM recusou o job $no (cota?): os ja submetidos seguem; submeta o resto depois" >&2
        break
    }
    echo "job ${j%%;*} ($((no + 1)) de $NOS)"
    echo "$(date '+%F %T') ${j%%;*} PARALELO=$PARALELO CASOS=$CASOS IREF=$IREF BLOCO=$BLOCO${ULTIMA:+ ULTIMA=$ULTIMA} SAIDA=$SAIDA" \
        >> logs/jobs.txt
done
echo "resultados em $SAIDA; acompanhe com ./status.sh"
