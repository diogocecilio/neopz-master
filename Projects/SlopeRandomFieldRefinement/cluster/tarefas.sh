#!/bin/bash
# tarefas.sh <rotulo>: um processo da fila de tarefas do Monte Carlo (o job.sbatch roda PARALELO deles por no).
# Tarefa = (caso, bloco): as realizacoes a..z = b BLOCO .. min((b + 1) BLOCO, N) - 1 do caso, numa rodada do EXE
# (SlopeRainMonteCarlo, SlopeDrawdownMonteCarlo ou outro com as opcoes n, primeira e saida) na pasta $SAIDA/<caso>
# (arquivos mc_<a>_<z>.csv, mc_mapa_<a>_<z>.csv e mc_<a>_<z>.log). Ordem: o bloco 0 de todos os casos, depois o 1, ...
# (os casos andam juntos). Cada processo pega a proxima tarefa livre criando o arquivo $SAIDA/<caso>/.tarefas/<a>_<z> com
# o noclobber do bash (open com O_EXCL: so um processo consegue, tambem entre nos; outro job pode entrar na fila a
# qualquer momento). Nao usa mkdir: o mkdir do coreutils em Rust (uutils, Ubuntu recente) nao e atomico, e dois processos
# pegavam a mesma tarefa. A tarefa esta feita quando existe mc_mapa_<a>_<z>.csv, gravado no fim da rodada. Se a rodada
# termina com erro, a trava ganha uma linha ERRO e a tarefa fica parada ate a limpeza do submeter.sh (sem jobs na fila).
# Nas opcoes dos casos, IREF vira $IREF e @SAIDA@ vira $SAIDA. A opcao prog=<nome> troca o programa do caso pelo
# <nome> da mesma pasta do EXE (ex.: prog=SlopeRainMonteCarlo numa campanha do Cho2010) e nao e passada a ele.
# Dependencia (amostragem por importancia): um caso com piloto=@SAIDA@/<id> e piloto_n=<n> nas opcoes so comeca quando
# os blocos das realizacoes 0..n-1 do caso <id> acabaram (mapa ou ERRO). Enquanto isso, o processo faz as outras tarefas;
# sem nenhuma livre, espera 60 s e olha de novo, ate ESPERA_MAX s (padrao 3 h) seguidos. O piloto de um caso escolhido
# entra na fila mesmo fora de CASOS, com o N do arquivo (sem ULTIMA): um bloco dele que precise ser refeito (ERRO limpo
# pelo submeter.sh) roda, e a amostragem nao fica esperando.
# Sem tarefa livre, ou com o arquivo $SAIDA/PARAR, o processo termina.
# Ambiente (job.sbatch): SAIDA, EXE, CASOS_ARQ, CASOS (todos, ids ou fases, separados por virgula), IREF, BLOCO, ULTIMA
# (opcional: N = ULTIMA + 1 em todos os casos escolhidos) e ESPERA_MAX (opcional)
set -uo pipefail
rotulo="$1"
LOG="$SAIDA/logs/p$rotulo.txt"
# todos os casos do arquivo; os escolhidos (CASOS) e os pilotos deles
todos_id=() todos_fase=() todos_n=() todos_opc=()
while read -r id fase n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    todos_id+=("$id") todos_fase+=("$fase") todos_n+=("$n") todos_opc+=("${opc//IREF/$IREF}")
done < "$CASOS_ARQ"
declare -A ESCOLHIDO PILOTO_EXTRA NCASO PRONTO
for i in "${!todos_id[@]}"; do
    if [ "$CASOS" = todos ]; then
        ESCOLHIDO[${todos_id[$i]}]=1
    else
        case ",$CASOS," in *",${todos_id[$i]},"* | *",fase${todos_fase[$i]},"*) ESCOLHIDO[${todos_id[$i]}]=1 ;; esac
    fi
done
for i in "${!todos_id[@]}"; do
    [ -n "${ESCOLHIDO[${todos_id[$i]}]:-}" ] || continue
    [[ " ${todos_opc[$i]} " =~ \ piloto=@SAIDA@/([^ ]+)\  ]] || continue
    [ -n "${ESCOLHIDO[${BASH_REMATCH[1]}]:-}" ] || PILOTO_EXTRA[${BASH_REMATCH[1]}]=1
done
ids=() opcs=() ns=() pil=() pnp=() exes=()
for i in "${!todos_id[@]}"; do
    id=${todos_id[$i]} n=${todos_n[$i]} opc=${todos_opc[$i]}
    NCASO[$id]=$n
    if [ -n "${ESCOLHIDO[$id]:-}" ]; then
        if [ -n "${ULTIMA:-}" ]; then n=$((ULTIMA + 1)) NCASO[$id]=$((ULTIMA + 1)); fi
    elif [ -z "${PILOTO_EXTRA[$id]:-}" ]; then
        continue
    fi
    p="" np=1000
    [[ " $opc " =~ \ piloto=@SAIDA@/([^ ]+)\  ]] && p="${BASH_REMATCH[1]}"
    [[ " $opc " =~ \ piloto_n=([0-9]+)\  ]] && np="${BASH_REMATCH[1]}"
    exe="$EXE"
    if [[ " $opc " =~ \ prog=([^ ]+)\  ]]; then
        exe="$(dirname "$EXE")/${BASH_REMATCH[1]}"
        opc="$(echo " $opc " | sed -e "s| prog=${BASH_REMATCH[1]} | |" -e 's/^ *//' -e 's/ *$//')"
    fi
    ids+=("$id") opcs+=("${opc//@SAIDA@/$SAIDA}") ns+=("$n") pil+=("$p") pnp+=("$np") exes+=("$exe")
done
for p in "${pil[@]}"; do
    [ -z "$p" ] || [ -n "${NCASO[$p]:-}" ] || { echo "$(date '+%F %T') piloto $p nao esta em $CASOS_ARQ: fim" >> "$LOG"; exit 1; }
done
echo "$(date '+%F %T') $(hostname -s) job ${SLURM_JOB_ID:-local}: ${#ids[@]} casos, blocos de $BLOCO" >> "$LOG"

# o piloto <id> (realizacoes 0..<n>-1) acabou? cada bloco dele tem o mapa ou ERRO na trava (uma vez pronto, fica pronto)
piloto_pronto() {
    [ -n "${PRONTO[$1]:-}" ] && return 0
    local d="$SAIDA/$1" np=$2 n=${NCASO[$1]} a z
    for ((a = 0; a < np && a < n; a += BLOCO)); do
        z=$(((a + BLOCO < n ? a + BLOCO : n) - 1))
        [ -f "$d/mc_mapa_${a}_${z}.csv" ] || grep -qs '^ERRO' "$d/.tarefas/${a}_${z}" || return 1
    done
    PRONTO[$1]=1
}

desde="" # inicio da espera pelos pilotos
while :; do
    espera=0 fez=0
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
            if [ -n "${pil[$i]}" ] && ! piloto_pronto "${pil[$i]}" "${pnp[$i]}"; then
                espera=1
                continue
            fi
            mkdir -p "$dir/.tarefas"
            # a trava: criar o arquivo, so se ainda nao existe (atomico); senao outro processo ja pegou
            ( set -o noclobber; echo "$(hostname -s) job ${SLURM_JOB_ID:-local} p$rotulo $(date '+%F %T')" \
                > "$dir/.tarefas/${a}_${z}" ) 2> /dev/null || continue
            fez=1 desde=""
            extra=()
            ((a == 0)) && extra=(det=1 figuras=-1)
            t0=$(date +%s)
            echo "$(date '+%F %T') ${ids[$i]} $a..$z: inicio" >> "$LOG"
            # shellcheck disable=SC2086  # as opcoes do caso se separam nos espacos
            (cd "$dir" && "${exes[$i]}" n=$((z - a + 1)) primeira=$a saida=mc ${opcs[$i]} ${extra[@]+"${extra[@]}"} \
                > "mc_${a}_${z}.log" 2>&1)
            c=$?
            if ((c == 0)); then
                echo "$(date '+%F %T') ${ids[$i]} $a..$z: fim em $(($(date +%s) - t0)) s" >> "$LOG"
            else
                echo "ERRO $(date '+%F %T') codigo $c" >> "$dir/.tarefas/${a}_${z}"
                echo "$(date '+%F %T') ${ids[$i]} $a..$z: ERRO (codigo $c), veja $dir/mc_${a}_${z}.log" >> "$LOG"
            fi
        done
        ((resta)) || break
    done
    ((espera)) || break # nada esperando um piloto: acabou
    ((fez)) && continue # fez alguma tarefa: olha a fila de novo
    if [ -z "$desde" ]; then
        desde=$(date +%s)
        echo "$(date '+%F %T') esperando os pilotos (as outras tarefas estao feitas ou com outros processos)" >> "$LOG"
    elif (($(date +%s) - desde > ${ESPERA_MAX:-10800})); then
        echo "$(date '+%F %T') os pilotos nao acabaram em ${ESPERA_MAX:-10800} s: fim" >> "$LOG"
        break
    fi
    sleep 60
done
echo "$(date '+%F %T') sem tarefas livres: fim" >> "$LOG"
