#!/bin/bash
# tarefas.sh <rotulo>: um processo da fila de tarefas do Monte Carlo da chuva (o job.sbatch roda PARALELO deles por no).
# Tarefa = (caso, bloco): as realizacoes a..z = b BLOCO .. min((b + 1) BLOCO, N) - 1 do caso, numa rodada do
# SlopeRainMonteCarlo na pasta $SAIDA/<caso> (arquivos mc_<a>_<z>.csv, mc_mapa_<a>_<z>.csv e mc_<a>_<z>.log). Ordem: o
# bloco 0 de todos os casos, depois o 1, ... (os casos andam juntos). Cada processo pega a proxima tarefa livre com
# mkdir $SAIDA/<caso>/.tarefas/<a>_<z>, que e atomico tambem entre nos (outro job pode entrar na fila a qualquer
# momento); a tarefa esta feita quando existe mc_mapa_<a>_<z>.csv, gravado no fim da rodada. Sem tarefa livre, ou com o
# arquivo $SAIDA/PARAR, o processo termina.
# Ambiente (job.sbatch): SAIDA, EXE, CASOS_ARQ, CASOS (todos, ids ou fases, separados por virgula), IREF, BLOCO e ULTIMA
# (opcional: N = ULTIMA + 1 em todos os casos escolhidos)
set -uo pipefail
rotulo="$1"
LOG="$SAIDA/logs/p$rotulo.txt"
ids=() opcs=() ns=()
while read -r id fase n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ]; then
        case ",$CASOS," in *",$id,"* | *",fase$fase,"*) ;; *) continue ;; esac
    fi
    [ -n "${ULTIMA:-}" ] && n=$((ULTIMA + 1))
    ids+=("$id") opcs+=("${opc//IREF/$IREF}") ns+=("$n")
done < "$CASOS_ARQ"
echo "$(date '+%F %T') $(hostname -s) job ${SLURM_JOB_ID:-local}: ${#ids[@]} casos, blocos de $BLOCO" >> "$LOG"

for ((b = 0; ; b++)); do
    resta=0
    for i in "${!ids[@]}"; do
        n=${ns[$i]} a=$((b * BLOCO))
        ((a < n)) || continue
        resta=1
        z=$(((a + BLOCO < n ? a + BLOCO : n) - 1))
        dir="$SAIDA/${ids[$i]}"
        [ -f "$dir/mc_mapa_${a}_${z}.csv" ] && continue
        if [ -e "$SAIDA/PARAR" ]; then
            echo "$(date '+%F %T') PARAR: fim" >> "$LOG"
            exit 0
        fi
        mkdir -p "$dir/.tarefas"
        mkdir "$dir/.tarefas/${a}_${z}" 2> /dev/null || continue # outro processo ja pegou
        echo "$(hostname -s) job ${SLURM_JOB_ID:-local} p$rotulo $(date '+%F %T')" > "$dir/.tarefas/${a}_${z}/dono"
        extra=()
        ((a == 0)) && extra=(det=1 figuras=-1)
        t0=$(date +%s)
        echo "$(date '+%F %T') ${ids[$i]} $a..$z: inicio" >> "$LOG"
        # shellcheck disable=SC2086  # as opcoes do caso se separam nos espacos
        if (cd "$dir" && "$EXE" n=$((z - a + 1)) primeira=$a saida=mc ${opcs[$i]} ${extra[@]+"${extra[@]}"} \
            > "mc_${a}_${z}.log" 2>&1); then
            echo "$(date '+%F %T') ${ids[$i]} $a..$z: fim em $(($(date +%s) - t0)) s" >> "$LOG"
        else
            echo "$(date '+%F %T') ${ids[$i]} $a..$z: ERRO (codigo $?), veja $dir/mc_${a}_${z}.log" >> "$LOG"
        fi
    done
    ((resta)) || break
done
echo "$(date '+%F %T') sem tarefas livres: fim" >> "$LOG"
