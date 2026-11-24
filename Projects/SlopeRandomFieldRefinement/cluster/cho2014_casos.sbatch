#!/bin/bash
# cho2014_casos.sbatch: a coluna 1D de talude infinito sob chuva de Cho (2014) (SlopeRainColumn, coluna.cpp; campo 1D
# lognormal de k_s pela base de Nystrom do seu lv, FS de Lu e Godt depois de cada passo, t_passagem localizado), com os
# casos do casos_cho2014.txt (os 15 da Fig. 10 e os estudos) e os FORM + line sampling do casos_cho2014_ls.txt (spec
# D.1 e D.2). Os jobs (submeter_cho2014.sh, que faz as bases de Nystrom no no de login antes) dividem o trabalho por
# travas em $SAIDA (noclobber, como no liu2024_casos.sbatch), em fases:
#   1. o deterministico de cada caso (k_s medio; ~10 s cada);
#   2. o Monte Carlo direto: min(N_mc do caso, NMC_MAX) realizacoes por caso em blocos de BLOCO_MC, intercalados entre
#      os casos (a realizacao n tem os mesmos xi em todos: uma semente so, a padrao do executavel);
#   3. para cada linha do casos_cho2014_ls.txt cujo caso ja tem o Monte Carlo pronto: o FORM por diferencas centrais
#      (form_kl.py prog=coluna gradiente=dif h=0.5 algoritmo=ihlrf, M termos, PARALELO analises em voo, ate FORM_IT
#      iteracoes), o line sampling na direcao dele (K = M termos, semente 7, ate CoV(P_f) <= COV, no minimo NMIN e no
#      maximo N_ls linhas) e a bissecao das 10 linhas de menor raiz (bissecao_ls.py). As tres etapas sao uma tarefa so
#      (uma trava); a retomada continua da etapa que falta (o line_sampling.py retoma do linhas.csv).
# Enquanto houver linha do LS esperando o MC do seu caso (com blocos rodando em outros jobs), o job espera em vez de sair.
# Retomada e parada (touch $SAIDA/PARAR: os jobs saem no fim da onda, o line sampling para de submeter linhas) como no
# liu2024_casos.sbatch; o submeter_cho2014.sh libera as travas de tarefas interrompidas (sem jobs cho14 na fila).
# Saidas: $SAIDA/<caso>/det/ (det_0_-1.csv, det.log), $SAIDA/<caso>/mc/ (mc_<a>_<b>.csv, log_<a>_<b>.txt, .tarefas/),
# $SAIDA/ls/<id>/ (form/form.json, form.log, ls/ls.json, ls.log, bis/bissecao.json, bis.log, FEITO_LS, FEITO) e
# $SAIDA/eventos.log. Resumo contra o artigo e a referencia em Python: resumo_cho2014.py.
# Uso: ./submeter_cho2014.sh [nos]. Ambiente: SAIDA ($DATA/cho2014), BASES ($SAIDA/bases), M (100), NY (80), EXTRA (opcoes
# do executavel a mais em todos os casos, ex. elementos=tri), PARALELO (24), COV (0.05), NMIN (30), N (600: o teto de
# linhas), FORM_IT (15), NMC_MAX (5000), BLOCO_MC (100), CASOS e LS_CASOS (todos, ou ids separados por virgulas; LS_CASOS=
# nenhum: sem FORM/LS), CASOS_TXT, LS_TXT, EXE (bin_coluna/SlopeRainColumn).
#SBATCH --job-name=cho14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
M="${M:-100}" NY="${NY:-80}"
SAIDA="${SAIDA:-$DADOS/cho2014}"
BASES="${BASES:-$SAIDA/bases}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}" FORM_IT="${FORM_IT:-15}"
NMC_MAX="${NMC_MAX:-5000}" BLOCO_MC="${BLOCO_MC:-100}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_cho2014.txt}" CASOS="${CASOS:-todos}"
LS_TXT="${LS_TXT:-$PWD/casos_cho2014_ls.txt}" LS_CASOS="${LS_CASOS:-todos}"
EXE="${EXE:-$PWD/bin_coluna/SlopeRainColumn}" PY="$(cd ../../Cho2010 && pwd)"
TEMPOS=6,12,16,18,20,22,24
NT=$(($(tr -cd , <<< "$TEMPOS" | wc -c) + 1)) # instantes por realizacao (linhas do csv)
# a base comum (spec D.1); o resto sao os padroes de producao do executavel (quad, massa consistente, newton 1e-6,80,1)
BASE="ny=$NY dt=0.1 dt_min=0.00625 M=$M tempos=$TEMPOS ${EXTRA:-}"
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_coluna.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }

# os casos: id, N do MC, lv e opcoes (CASOS: so os ids pedidos)
IDS=() OPCS=() NMCS=() LVS=()
while read -r id pf18 pf24 nmc opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ] && [[ ",$CASOS," != *",$id,"* ]]; then continue; fi
    lv="$(sed -n 's/.*\blv=\([0-9.]*\).*/\1/p' <<< "$opc")"
    [ -n "$lv" ] || { echo "caso $id sem lv= nas opcoes" >&2; exit 1; }
    [ "$nmc" -gt "$NMC_MAX" ] && nmc="$NMC_MAX"
    IDS+=("$id") OPCS+=("$opc") NMCS+=("$nmc") LVS+=("$lv")
done < "$CASOS_TXT"
[ ${#IDS[@]} -gt 0 ] || { echo "nenhum caso (CASOS=$CASOS em $CASOS_TXT)" >&2; exit 1; }
indice() { local i; for i in "${!IDS[@]}"; do [ "${IDS[$i]}" = "$1" ] && echo "$i" && return 0; done; return 1; }
# as linhas do FORM + LS: id, o indice do caso, T, falha e o teto de linhas (so as de casos presentes)
LIDS=() LCASO=() LT=() LFALHA=() LN=()
if [ "$LS_CASOS" != nenhum ]; then
    while read -r id caso t falha nls; do
        [[ -z "$id" || "$id" == \#* ]] && continue
        if [ "$LS_CASOS" != todos ] && [[ ",$LS_CASOS," != *",$id,"* ]]; then continue; fi
        i="$(indice "$caso")" || continue
        [ "$nls" -gt "$N" ] && nls="$N"
        LIDS+=("$id") LCASO+=("$i") LT+=("$t") LFALHA+=("$falha") LN+=("$nls")
    done < "$LS_TXT"
fi

# pega a trava <arquivo> (cria so se nao existe: noclobber, atomico tambem entre nos). 0: pegou; 1: outro pegou. Uma
# trava VAZIA com mais de 10 min e lixo do NFS ("Stale file handle"): sai e tenta de novo (como no tarefas.sh)
pega_trava() {
    local t
    for t in 1 2 3 4 5; do
        mkdir -p "$(dirname "$1")" 2> /dev/null
        ( set -o noclobber; echo "$(hostname -s) job $JOB $(date '+%F %T')" > "$1" ) 2> /dev/null && return 0
        if [ -e "$1" ]; then
            if [ -n "$(find "$1" -maxdepth 0 -empty -mmin +10 2> /dev/null)" ]; then
                rm -f "$1"
                continue
            fi
            return 1
        fi
        sleep 1
    done
    return 1
}
com_erro() { [ -e "$1" ] && grep -q '^ERRO' "$1"; }
base_de() { echo "$BASES/base_lv${LVS[$1]}_M$M.txt"; } # a base de Nystrom do caso <i>
det_feito() { [ -s "$SAIDA/$1/det/det_0_-1.csv" ] && [ "$(wc -l < "$SAIDA/$1/det/det_0_-1.csv")" -gt "$NT" ]; }
mc_feito() { # <pasta mc> <a> <b>: o csv do bloco com as (b - a + 1) x NT linhas
    local f="$1/mc_$2_$3.csv"
    [ -s "$f" ] && [ $(($(wc -l < "$f") - 1)) -ge $((($3 - $2 + 1) * NT)) ]
}
blocos() { # <i>: os blocos "<a> <b>" do caso
    local a b
    for ((a = 0; a < NMCS[$1]; a += BLOCO_MC)); do
        b=$((a + BLOCO_MC - 1)); [ "$b" -ge "${NMCS[$1]}" ] && b=$((NMCS[$1] - 1))
        echo "$a $b"
    done
}
# o MC do caso <i>: 0 pronto; 1 falta (algum bloco livre ou rodando); 2 algum bloco com ERRO (o LS dele nao roda)
mc_estado() {
    local a b d="$SAIDA/${IDS[$1]}/mc" r=0
    while read -r a b; do
        mc_feito "$d" "$a" "$b" && continue
        com_erro "$d/.tarefas/${a}_${b}" && return 2
        r=1
    done < <(blocos "$1")
    return "$r"
}

# a base de Nystrom do caso <i> no no: copia local (as analises do no leem a mesma), ou a de $BASES se a copia falhar.
# Feita uma vez no inicio do job (NYSL), antes de qualquer tarefa em paralelo
nys_local() {
    local b loc
    b="$(base_de "$1")"
    loc="${TMPDIR:-/tmp}/cho2014_$(basename "$b")"
    if [ ! -s "$loc" ] || [ "$(stat -c %s "$loc")" != "$(stat -c %s "$b")" ]; then
        cp "$b" "$loc.$JOB" 2> /dev/null && mv -f "$loc.$JOB" "$loc" || loc="$b"
    fi
    echo "$loc"
}
NYSL=()

# fase 1: os deterministicos (todos ao mesmo tempo, ate PARALELO)
fase_det() {
    local i id d n=0
    for i in "${!IDS[@]}"; do
        [ "$n" -ge "$PARALELO" ] && break
        id="${IDS[$i]}" d="$SAIDA/${IDS[$i]}"
        det_feito "$id" && continue
        pega_trava "$d/.det.trava" || continue
        (
            nys="${NYSL[$i]}"
            mkdir -p "$d/det" && cd "$d/det" || exit 1
            # shellcheck disable=SC2086
            "$EXE" n=0 det=1 $BASE ${OPCS[$i]} nystrom="$nys" saida=det > det.log 2>&1
            if det_feito "$id"; then
                diz "deterministico $id: FS $(awk -F, 'NR > 1 { printf "%s h %s; ", $2, $3 }' det_0_-1.csv)t_passagem" \
                    "$(awk -F, 'END { print $8 }' det_0_-1.csv) h"
            else
                echo "ERRO $(date '+%F %T')" >> "$d/.det.trava"
                diz "deterministico $id FALHOU (veja $d/det/det.log)"
            fi
        ) &
        n=$((n + 1))
    done
    wait
}

# fase 2: os blocos do Monte Carlo, intercalados entre os casos: "<i> <a> <b>"
mc_tarefas() {
    local i a b maior=0
    for i in "${!IDS[@]}"; do [ "${NMCS[$i]}" -gt "$maior" ] && maior="${NMCS[$i]}"; done
    for ((a = 0; a < maior; a += BLOCO_MC)); do
        for i in "${!IDS[@]}"; do
            [ "$a" -lt "${NMCS[$i]}" ] || continue
            b=$((a + BLOCO_MC - 1)); [ "$b" -ge "${NMCS[$i]}" ] && b=$((NMCS[$i] - 1))
            echo "$i $a $b"
        done
    done
}
# pega ate PARALELO blocos livres, roda todos e espera; 1 se nao havia bloco livre
mc_onda() {
    local n=0 i a b d nys
    while read -r i a b; do
        [ "$n" -ge "$PARALELO" ] && break
        d="$SAIDA/${IDS[$i]}/mc"
        mc_feito "$d" "$a" "$b" && continue
        pega_trava "$d/.tarefas/${a}_${b}" || continue
        rm -f "$d/mc_${a}_${b}.csv"
        nys="${NYSL[$i]}"
        (
            cd "$d" || exit 1
            # shellcheck disable=SC2086
            "$EXE" n=$((b - a + 1)) primeira="$a" $BASE ${OPCS[$i]} nystrom="$nys" saida=mc > "log_${a}_${b}.txt" 2>&1 \
                && mc_feito "$d" "$a" "$b" || echo "ERRO $(date '+%F %T')" >> ".tarefas/${a}_${b}"
        ) &
        n=$((n + 1))
    done < <(mc_tarefas)
    [ "$n" -gt 0 ] || return 1
    wait
    return 0
}

# fase 3: o FORM, o line sampling e a bissecao da linha <j> do LS (com a trava ja pega)
ls_item() {
    local j="$1" id="${LIDS[$1]}" i="${LCASO[$1]}" d="$SAIDA/ls/${LIDS[$1]}" nys pid rc
    nys="${NYSL[$i]}"
    mkdir -p "$d"
    if [ ! -s "$d/form/form.json" ]; then
        diz "FORM $id: caso ${IDS[$i]} (${OPCS[$i]}), T ${LT[$j]} h, falha ${LFALHA[$j]}"
        # shellcheck disable=SC2086
        if python3 -I "$PY/form_kl.py" "$EXE" "$d/form" prog=coluna falha="${LFALHA[$j]}" T="${LT[$j]}" gradiente=dif \
            h=0.5 algoritmo=ihlrf iteracoes="$FORM_IT" paralelos="$PARALELO" $BASE ${OPCS[$i]} nystrom="$nys" \
            > "$d/form.log" 2>&1; then
            diz "FORM $id pronto: $(grep -m 1 'beta = ' "$d/form.log" | tr -s ' ')$(grep -q 'NAO convergiu' "$d/form.log" && echo ' (nao convergiu)')"
        else
            echo "ERRO form $(date '+%F %T')" >> "$d/.trava"
            diz "FORM $id FALHOU (veja $d/form.log)"
            return
        fi
    fi
    [ -e "$SAIDA/PARAR" ] && { rm -f "$d/.trava"; return; }
    if [ ! -e "$d/FEITO_LS" ]; then
        rm -f "$d/ls/PARAR"
        diz "line sampling $id"
        python3 -I "$PY/line_sampling.py" "$EXE" "$d/form/form.json" "$d/ls" K="$M" semente=7 cov="$COV" nmin="$NMIN" \
            n="${LN[$j]}" paralelos="$PARALELO" threads=0 nystrom="$nys" >> "$d/ls.log" 2>&1 &
        pid=$!
        while kill -0 "$pid" 2> /dev/null; do
            [ -e "$SAIDA/PARAR" ] && [ ! -e "$d/ls/PARAR" ] && touch "$d/ls/PARAR"
            sleep 30
        done
        wait "$pid"
        rc=$?
        if [ "$rc" -ne 0 ]; then
            echo "ERRO ls $(date '+%F %T')" >> "$d/.trava"
            diz "line sampling $id FALHOU (codigo $rc; veja $d/ls.log)"
            return
        fi
        if [ -s "$d/ls/ls.json" ] && ! grep -q '"motivo": "PARAR' "$d/ls/ls.json"; then
            touch "$d/FEITO_LS"
            diz "line sampling $id pronto: $(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1]))
print("P_f %.4e, CoV %.1f%%, %d linhas (%d perdidas), %d analises" % (d["P_f"], 100*(d["CoV"] or 0), d["linhas"], d["perdidas"], d["analises"]) if d["P_f"] is not None else "sem P_f")' "$d/ls/ls.json")"
        else
            diz "line sampling $id parado (PARAR)"
            rm -f "$d/.trava" # outro job retoma depois
            return
        fi
    fi
    if [ ! -s "$d/bis/bissecao.json" ]; then
        diz "bissecao $id"
        if python3 -I "$PY/bissecao_ls.py" "$EXE" "$d/ls" "$d/bis" linhas=10 paralelos="$PARALELO" > "$d/bis.log" 2>&1; then
            diz "bissecao $id pronta: $(tail -n 1 "$d/bis.log")"
        else
            echo "ERRO bis $(date '+%F %T')" >> "$d/.trava"
            diz "bissecao $id FALHOU (veja $d/bis.log)"
            return
        fi
    fi
    touch "$d/FEITO"
    rm -f "$d/.trava"
}

# a primeira linha do LS livre (o MC do caso pronto, sem FEITO e sem trava)
livre_ls() {
    local j
    for j in "${!LIDS[@]}"; do
        [ -e "$SAIDA/ls/${LIDS[$j]}/FEITO" ] && continue
        [ -e "$SAIDA/ls/${LIDS[$j]}/.trava" ] && continue
        mc_estado "${LCASO[$j]}" || continue
        echo "$j"
        return 0
    done
    return 1
}
# alguma linha do LS esperando o MC do seu caso (sem trava, com o MC sem ERRO e ainda por terminar)?
espera_ls() {
    local j r
    for j in "${!LIDS[@]}"; do
        [ -e "$SAIDA/ls/${LIDS[$j]}/FEITO" ] && continue
        [ -e "$SAIDA/ls/${LIDS[$j]}/.trava" ] && continue
        mc_estado "${LCASO[$j]}"
        r=$?
        [ "$r" = 1 ] && return 0
    done
    return 1
}

diz "inicio em $(hostname -s): ${#IDS[@]} casos, ${#LIDS[@]} FORM/LS, M $M, ny $NY, PARALELO $PARALELO, MC ate" \
    "$NMC_MAX por caso (blocos de $BLOCO_MC), CoV alvo $COV, linhas $NMIN a $N, FORM ate $FORM_IT iteracoes," \
    "executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
for i in "${!IDS[@]}"; do
    [ -s "$(base_de "$i")" ] || { diz "falta a base $(base_de "$i") (o submeter_cho2014.sh faz): saindo"; exit 1; }
    NYSL[$i]="$(nys_local "$i")"
done
while :; do
    [ -e "$SAIDA/PARAR" ] && { diz "achei $SAIDA/PARAR: saindo"; break; }
    fase_det
    if j="$(livre_ls)"; then
        pega_trava "$SAIDA/ls/${LIDS[$j]}/.trava" || continue
        ls_item "$j"
        continue
    fi
    mc_onda && continue
    if espera_ls; then # o MC de um caso do LS ainda roda em outro job: espera para fazer o LS
        sleep 60
        continue
    fi
    diz "nada mais livre: saindo"
    break
done
diz "fim"
