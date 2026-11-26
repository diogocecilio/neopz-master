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
# Nas opcoes dos casos, IREF vira $IREF e @SAIDA@ vira $SAIDA. Opcoes da fila, que nao vao para o programa:
#   prog=<nome>      o programa do caso e o <nome> da mesma pasta do EXE (ex.: prog=SlopeRainMonteCarlo numa campanha
#                    do Cho2010)
#   espera=<a>,<b>   o caso so comeca quando os arquivos existem (ex.: o ponto de projeto que o form.sh grava); com um
#                    <arquivo>.ERRO (quem o fazia desistiu), o caso fica de fora
#   ref=<id>, par=<id>  o Monte Carlo direto de referencia e o caso da comparacao pareada (so para o baixar.sh)
#   cov_alvo=<cv>    parada pelo CoV (amostragem por importancia ou MC direto, log_peso = 0; com varios instantes, o menor
#                    FS de cada realizacao; com crpc= (uma intensidade), a falha da primeira passagem do mc_crpc_*.csv):
#                    antes de pegar um bloco novo do caso, a P_f ponderada
#                    sum(I w) / k e o seu CoV sao calculados nas realizacoes 0..k-1 dos blocos ja prontos sem buraco (na
#                    ordem; I = FS_min < 1, w = exp(log_peso), so as linhas com situacao ok e FS_min > 0, como o
#                    amostragem_importancia.py; no GIM, Gamma = 0 conta como falha). Com CoV <= cv, k >= cov_n (500) e pelo menos cov_falhas (100) falhas
#                    sorteadas, o caso grava <caso>/ATINGIU (k, falhas, P_f, CoV, data) e nao lanca mais blocos; os que ja
#                    rodam terminam. O N do caso vira so o teto. O N em que o CoV cruza o alvo sai depois, das realizacoes
#                    na ordem (Projects/Cho2010/ismc_parada.py)
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
ids=() opcs=() ns=() pil=() pnp=() exes=() esps=() cvas=() cvns=() cvfs=()
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
    # as opcoes da fila saem daqui
    exe="$EXE" esp="" resto="" cva="" cvn=500 cvf=100
    set -f
    for w in $opc; do
        case "$w" in
        prog=*) exe="$(dirname "$EXE")/${w#prog=}" ;;
        espera=*) esp="${w#espera=}" ;;
        cov_alvo=*) cva="${w#cov_alvo=}" ;;
        cov_n=*) cvn="${w#cov_n=}" ;;
        cov_falhas=*) cvf="${w#cov_falhas=}" ;;
        ref=* | par=*) ;;
        *) resto+="${resto:+ }$w" ;;
        esac
    done
    set +f
    ids+=("$id") opcs+=("${resto//@SAIDA@/$SAIDA}") ns+=("$n") pil+=("$p") pnp+=("$np") exes+=("$exe")
    esps+=("${esp//@SAIDA@/$SAIDA}") cvas+=("$cva") cvns+=("$cvn") cvfs+=("$cvf")
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

# os arquivos de espera= (separados por virgula): 0 se existem todos, 1 se falta algum, 2 se algum tem o <arquivo>.ERRO
arquivos_prontos() {
    local f r=0
    local -a fs
    IFS=, read -ra fs <<< "$1"
    for f in "${fs[@]}"; do
        [ -e "$f.ERRO" ] && return 2
        [ -e "$f" ] || r=1
    done
    return $r
}
declare -A DE_FORA

# o caso i ja chegou ao CoV alvo (cov_alvo=)? 0: sim (ATINGIU gravado agora ou antes); 1: ainda nao. As somas do
# prefixo pronto (k, falhas, sum I w, sum (I w)^2) ficam guardadas por processo: cada bloco e lido uma vez so
declare -A ATINGIU_CASO PREF_A PREF_S
atingiu() {
    local i=$1 dir="$SAIDA/${ids[$1]}" n=${ns[$1]} a z linha k nf s1 s2 pf cv ok
    local -a novos=()
    [ -n "${ATINGIU_CASO[$i]:-}" ] && return 0
    [ -f "$dir/ATINGIU" ] && { ATINGIU_CASO[$i]=1; return 0; }
    for ((a = ${PREF_A[$i]:-0}; a < n; a += BLOCO)); do
        z=$(((a + BLOCO < n ? a + BLOCO : n) - 1))
        [ -f "$dir/mc_mapa_${a}_${z}.csv" ] || break
        if [ -f "$dir/mc_crpc_${a}_${z}.csv" ]; then novos+=("$dir/mc_crpc_${a}_${z}.csv") # primeira passagem (crpc=)
        elif [ -f "$dir/mc_${a}_${z}.csv" ]; then novos+=("$dir/mc_${a}_${z}.csv"); fi
    done
    ((${#novos[@]})) || return 1 # nenhum bloco novo no prefixo desde a ultima conta
    PREF_A[$i]=$a
    linha="$(awk -F, -v ini="${PREF_S[$i]:-0 0 0 0}" -v alvo="${cvas[$i]}" -v nmin="${cvns[$i]}" -v fmin="${cvfs[$i]}" '
        BEGIN { split(ini, s, " "); k = s[1] + 0; nf = s[2] + 0; s1 = s[3] + 0; s2 = s[4] + 0 }
        FNR == 1 { crp = 0
                   for (c = 1; c <= NF; c++) { if ($c == "falha") cfa = c; if ($c == "DC_h") crp = 1
                                                 if ($c == "real") cr = c; if ($c == "FS_min") cf = c
                                                 if ($c == "log_peso") cw = c; if ($c == "situacao") cs = c
                                                 if ($c == "metodo") cm = c }
                   next }
        # primeira passagem (crpc= com uma intensidade): falha = "sim"; "nao" sobrevive; descartada e erros ficam fora
        crp { if ($cr + 0 >= 0 && ($cfa == "sim" || $cfa == "nao")) {
                  r = FILENAME SUBSEP $cr; vv = ($cfa == "sim") ? 0.5 : 2
                  if (!(r in fs) || vv < fs[r]) fs[r] = vv
                  lw[r] = 0 }
              next }
        # GIM com Gamma = 0 (situacao fator, sem equilibrio nem com 1% da gravidade): falha de verdade, entra
        { ok = ($cs == "ok" && $cf + 0 > 0) || (cm && $cm == "gim" && $cs == "fator" && $cf ~ /^0(\.0*)?$/) }
        # uma realizacao por (arquivo, real): com varios instantes (tempos=24,30), o menor FS das linhas ok dela
        $cr + 0 >= 0 && ok {
            r = FILENAME SUBSEP $cr
            if (!(r in fs) || $cf + 0 < fs[r]) fs[r] = $cf + 0
            lw[r] = (cw ? $cw + 0 : 0) }
        END { for (r in fs) { k++; if (fs[r] < 1) { iw = exp(lw[r]); s1 += iw; s2 += iw * iw; nf++ } }
              pf = (k > 0 ? s1 / k : 0); cv = 9
              if (k > 1 && pf > 0) { v = (s2 - k * pf * pf) / (k - 1); cv = (v > 0 ? sqrt(v / k) : 0) / pf }
              printf "%d %d %.17g %.17g %.6g %.4f %d\n", k, nf, s1, s2, pf, cv, (k >= nmin && nf >= fmin && cv <= alvo) }' \
        "${novos[@]}")" || return 1
    read -r k nf s1 s2 pf cv ok <<< "$linha"
    PREF_S[$i]="$k $nf $s1 $s2"
    ((ok)) || return 1
    echo "$k $nf $pf $cv $(date '+%F %T')" > "$dir/ATINGIU"
    echo "$(date '+%F %T') ${ids[$i]}: CoV alvo ${cvas[$i]} atingido com k = $k (falhas $nf, P_f $pf, CoV $cv)" >> "$LOG"
    ATINGIU_CASO[$i]=1
    return 0
}

# pega a trava <arquivo>: cria o arquivo so se ele ainda nao existe (noclobber, atomico). 0: pegou; 1: outro processo
# pegou. Se a criacao falha e o arquivo nao existe (o sistema de arquivos falhou, ou a pasta .tarefas ainda estava sendo
# criada por outro processo), tenta de novo, ate 5 vezes: antes, essa falha passava por "outro pegou", e o bloco ficava
# sem ninguem ate o fim da passada (na etapa1, o bloco 10..19 do piloto do r05_ai)
pega_trava() {
    local t
    for t in 1 2 3 4 5; do
        mkdir -p "$(dirname "$1")" 2> /dev/null
        ( set -o noclobber; echo "$(hostname -s) job ${SLURM_JOB_ID:-local} p$rotulo $(date '+%F %T')" > "$1" ) 2> /dev/null &&
            return 0
        if [ -e "$1" ]; then
            # uma trava VAZIA e com mais de 10 min nunca foi de ninguem: o sistema de arquivos de rede criou o arquivo e
            # devolveu erro (na etapa1, "Stale file handle" no bloco 10..19, parado 7 h). Remove e tenta de novo
            if [ -n "$(find "$1" -maxdepth 0 -empty -mmin +10 2> /dev/null)" ]; then
                rm -f "$1"
                echo "$(date '+%F %T') trava vazia e antiga removida: $1" >> "$LOG"
                continue
            fi
            return 1
        fi
        echo "$(date '+%F %T') a trava $1 nao foi criada (tentativa $t)" >> "$LOG"
        sleep 1
    done
    return 1
}

desde="" # inicio da espera pelos pilotos e arquivos
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
            if [ -n "${esps[$i]}" ]; then
                arquivos_prontos "${esps[$i]}"
                case $? in
                1) espera=1; continue ;;
                2) [ -n "${DE_FORA[$i]:-}" ] || echo "$(date '+%F %T') ${ids[$i]}: ha um .ERRO em ${esps[$i]}: fica de fora" >> "$LOG"
                    DE_FORA[$i]=1
                    continue ;;
                esac
            fi
            # parada pelo CoV: o caso que ja chegou ao alvo nao lanca mais blocos
            [ -n "${cvas[$i]}" ] && atingiu "$i" && continue
            pega_trava "$dir/.tarefas/${a}_${z}" || continue
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
            break 2 # depois de cada tarefa, a busca recomeca do bloco 0: um bloco pulado e os pilotos vem primeiro
        done
        ((resta)) || break
    done
    ((fez)) && continue # fez uma tarefa: olha a fila de novo, do bloco 0
    ((espera)) || break # nada feito e nada esperando um piloto ou arquivo: acabou
    if [ -z "$desde" ]; then
        desde=$(date +%s)
        echo "$(date '+%F %T') esperando os pilotos ou os arquivos de espera= (as outras tarefas estao feitas ou com outros" \
            "processos)" >> "$LOG"
    elif (($(date +%s) - desde > ${ESPERA_MAX:-10800})); then
        echo "$(date '+%F %T') os pilotos ou arquivos nao ficaram prontos em ${ESPERA_MAX:-10800} s: fim" >> "$LOG"
        break
    fi
    sleep 60
done
echo "$(date '+%F %T') sem tarefas livres: fim" >> "$LOG"
