#!/bin/bash
# hu2024_casos.sbatch: o caso A, Hu et al. (2024) (talude de 10 m a 45 graus sob chuva, c' e phi' lognormais com
# r(c, phi) = -0,5, k_s fixo), com o SlopeRainMonteCarlo: FORM com o gradiente modal, line sampling e Monte Carlo direto
# por cenario de chuva (casos_hu2024.txt), como o ls_ceron_casos.sbatch. Com k_s fixo, um transiente por cenario serve a
# todas as realizacoes: ele e calculado uma vez (so_transiente=1 grava_pressoes=) e todas as analises o leem (pressoes=).
# Os jobs (submeter_hu2024.sh) dividem o trabalho por travas em $SAIDA (noclobber, como no tarefas.sh), em fases:
#   0. a base de Nystrom do talude (M termos, nucleo exp 20/2 m, grade 0,3 x 0,1 m; nystrom_base.py com o python do
#      venv) e o transiente de cada caso (malha TRANS, ate 72 h, p_w de 6 em 6 h): ate PARALELO tarefas por job;
#   1. o deterministico (as medias, FS de 6 em 6 h de 0 a 72 h) e o FORM (M termos por campo, gradiente modal, ate
#      FORM_IT iteracoes de 1 analise com FORM_THREADS threads) de cada caso;
#   2. o line sampling de cada caso com o FORM pronto, um caso por vez com PARALELO linhas em voo, ate CoV(P_f) <= COV
#      (no minimo NMIN e no maximo o teto do caso, limitado a N linhas), semente 7 em todos os casos;
#   3. o Monte Carlo direto: min(N_mc do caso, NMC_MAX) realizacoes por caso em blocos de BLOCO_MC, intercalados entre os
#      casos (a realizacao n tem os mesmos numeros aleatorios em todos os casos), e a conferencia da primeira passagem:
#      as NCHECK primeiras realizacoes de cada caso com o FS em todos os instantes (mc13/).
# Instantes do FS (TEMPOS): fim (padrao) = o fim da chuva D (com refinamento adaptativo) e D + 6 h (na malha de D): o
# deterministico e o piloto de 22/11 dao o menor FS em t = D (o artigo: "the FOS reaches minimum value when the rainfall
# stop"), e a conferencia mc13/ mede isso em cada caso; todos = de 6 em 6 h de 0 a 72 h, todos adaptativos (13 SRM por
# analise). A falha e o menor FS dos instantes abaixo de 1 (primeira passagem).
# Retomada: o FORM interrompido recomeca do zero; o line sampling continua do linhas.csv; o bloco do MC incompleto e
# refeito. As travas de um job que morreu sao liberadas pelo submeter_hu2024.sh quando nao ha jobs huA na fila.
# Parar sem perder nada: touch $SAIDA/PARAR (os line samplings param de submeter linhas e terminam as em voo; os blocos
# do MC em curso terminam). Para voltar: rm $SAIDA/PARAR e ./submeter_hu2024.sh de novo.
# Saidas por caso em $SAIDA/<id>/: pressoes.txt (o transiente), det/ (det_0_-1.csv), form/ (form.json), form.log, ls/
# (ls.json, linhas.csv, analises.csv, perdidas.csv), ls.log, FEITO_LS, mc/ (mc_<a>_<b>.csv), mc13/ e as travas. Eventos
# em $SAIDA/eventos.log. Resumo contra o artigo: resumo_hu2024.py.
# Uso: ./submeter_hu2024.sh [nos]. Ambiente: KS (2e-5: o k_s do texto do artigo; 5.2e-6: o patamar da Fig. 4b), SAIDA
# ($DATA/hu2024_ks$KS), M (150), BASE_NYS ($DATA/hu2024/base_hu_M$M.txt), TRANS (tri4), TEMPOS (fim), PARALELO (24), COV
# (0.05), NMIN (30), N (600), FORM_IT (10), FORM_THREADS (4), NMC_MAX (8200), BLOCO_MC (10), NCHECK (20), CASOS (todos, ou
# ids separados por virgulas), EXE (bin_hu/SlopeRainMonteCarlo), PYV ($HOME/venv_la/bin/python, com numpy e scipy),
# TODOS (os instantes do transiente, do deterministico e da conferencia: de 6 em 6 h de 0 a 72 h).
#SBATCH --job-name=huA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
KS="${KS:-2e-5}" M="${M:-150}"
SAIDA="${SAIDA:-$DADOS/hu2024_ks$KS}"
BASE_NYS="${BASE_NYS:-$DADOS/hu2024/base_hu_M$M.txt}"
TRANS="${TRANS:-tri4}" TEMPOS="${TEMPOS:-fim}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}" FORM_IT="${FORM_IT:-10}" FORM_THREADS="${FORM_THREADS:-4}"
NMC_MAX="${NMC_MAX:-8200}" BLOCO_MC="${BLOCO_MC:-10}" NCHECK="${NCHECK:-20}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_hu2024.txt}" CASOS="${CASOS:-todos}"
EXE="${EXE:-$PWD/bin_hu/SlopeRainMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)" PYV="${PYV:-$HOME/venv_la/bin/python}"
TALUDE=0,0,70,0,70,10,40,10,30,20,0,20
# o caso A: geometria, nivel 5 m abaixo do pe e horizontal, succao ate 50 kPa, lateral esquerda fechada e direita com
# carga (5 m), van Genuchten-Mualem da Fig. 4, chi = S_e, gamma, c', phi' e os campos (Tabela 1), malhas
BASE="geometria=$TALUDE hw=5 na_beta=0 succao_max=50 esquerda=fechada direita=carga vg=0.399,0.024,0.0357,1.71 chi=se"
BASE+=" gamma=17.8 c=8 phi=24 covc=0.15 covphi=0.1 covk=0 rcphi=-0.5 ks=$KS transiente=$TRANS fator=tri1 nref=4"
TODOS="${TODOS:-0,6,12,18,24,30,36,42,48,54,60,66,72}" # os instantes do transiente, do deterministico e do mc13
NTODOS=$(($(tr -cd , <<< "$TODOS" | wc -c) + 1))
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_hu2024.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }

# os casos: id, P_f do artigo, N do MC, teto de linhas do LS e opcoes (CASOS: so os ids pedidos)
IDS=() OPCS=() NMCS=() NLS=() DS=()
while read -r id pf nmc nls opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ] && [[ ",$CASOS," != *",$id,"* ]]; then continue; fi
    d="$(sed -n 's/.*chuva=\([0-9.]*\).*/\1/p' <<< "$opc")"
    [ -n "$d" ] || { echo "caso $id sem chuva= nas opcoes" >&2; exit 1; }
    [ "$nmc" -gt "$NMC_MAX" ] && nmc="$NMC_MAX"
    [ "$nls" -gt "$N" ] && nls="$N"
    IDS+=("$id") OPCS+=("$opc") NMCS+=("$nmc") NLS+=("$nls") DS+=("$d")
done < "$CASOS_TXT"
[ ${#IDS[@]} -gt 0 ] || { echo "nenhum caso (CASOS=$CASOS em $CASOS_TXT)" >&2; exit 1; }
# os instantes do caso <i> (tempos= e adapta=)
instantes() {
    local d="${DS[$1]}"
    if [ "$TEMPOS" = todos ]; then echo "tempos=$TODOS"; else echo "tempos=$d,$(awk -v d="$d" 'BEGIN { print d + 6 }') adapta=$d"; fi
}

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
base_feita() { [ -s "$BASE_NYS" ]; }
tr_feito() { [ -s "$SAIDA/$1/pressoes.txt" ]; }
form_feito() { [ -s "$SAIDA/$1/form/form.json" ]; }
det_feito() { [ -s "$SAIDA/$1/det/det_0_-1.csv" ] && [ "$(wc -l < "$SAIDA/$1/det/det_0_-1.csv")" -gt "$NTODOS" ]; }
ls_feito() { [ -e "$SAIDA/$1/FEITO_LS" ]; }
mc_feito() { # <pasta mc> <a> <b> <instantes por realizacao>: o csv do bloco com as (b - a + 1) x instantes linhas
    local f="$1/mc_$2_$3.csv"
    [ -s "$f" ] && [ $(($(wc -l < "$f") - 1)) -ge $((($3 - $2 + 1) * $4)) ]
}

# a base de Nystrom no no: copia local (todas as analises do no leem a mesma), ou a de $DATA se a copia falhar
NYS_LOCAL=""
base_local() {
    local loc="${TMPDIR:-/tmp}/hu2024_$(basename "$BASE_NYS")"
    if [ ! -s "$loc" ] || [ "$(stat -c %s "$loc")" != "$(stat -c %s "$BASE_NYS")" ]; then
        cp "$BASE_NYS" "$loc.$JOB" 2> /dev/null && mv -f "$loc.$JOB" "$loc" || loc="$BASE_NYS"
    fi
    NYS_LOCAL="$loc"
}

# fase 0: a base (uma tarefa so) e os transientes; pega ate PARALELO tarefas livres, roda todas e espera
fase0() {
    local i id d n=0
    if ! base_feita && pega_trava "$BASE_NYS.trava"; then
        diz "base de Nystrom $BASE_NYS (M = $M)"
        (
            mkdir -p "$(dirname "$BASE_NYS")"
            if "$PYV" -I "$PY/nystrom_base.py" "$BASE_NYS.tmp" M="$M" talude="$TALUDE" dx=0.3 dy=0.1 \
                > "$BASE_NYS.log" 2>&1 && mv -f "$BASE_NYS.tmp" "$BASE_NYS"; then
                diz "base pronta: $(grep -m 1 "nos de" "$BASE_NYS.log" | cut -d, -f1); $(grep -m 1 "M = *$M:" "$BASE_NYS.log" | tr -s ' ')"
            else
                echo "ERRO $(date '+%F %T')" >> "$BASE_NYS.trava"
                diz "base de Nystrom FALHOU (veja $BASE_NYS.log)"
            fi
        ) &
        n=$((n + 1))
    fi
    for i in "${!IDS[@]}"; do
        [ "$n" -ge "$PARALELO" ] && break
        id="${IDS[$i]}" d="$SAIDA/${IDS[$i]}"
        tr_feito "$id" && continue
        pega_trava "$d/.tr.trava" || continue
        diz "transiente $id (${OPCS[$i]})"
        (
            # shellcheck disable=SC2086
            if PZ_THREADS=4 "$EXE" $BASE ${OPCS[$i]} tempos=$TODOS so_transiente=1 grava_pressoes="$d/pressoes.txt" \
                saida="$d/tr" > "$d/transiente.log" 2>&1; then
                diz "transiente $id pronto: $(grep -m 1 'transiente ate' "$d/transiente.log" | tr -s ' ')"
            else
                echo "ERRO $(date '+%F %T')" >> "$d/.tr.trava"
                diz "transiente $id FALHOU (veja $d/transiente.log)"
            fi
        ) &
        n=$((n + 4))
    done
    wait
}

# fase 1: deterministico e FORM (com o transiente e a base prontos)
fase1() {
    local tipo i id d n=0
    base_feita || return 0
    base_local
    for tipo in form det; do
        for i in "${!IDS[@]}"; do
            [ "$n" -ge "$PARALELO" ] && break 2
            id="${IDS[$i]}" d="$SAIDA/${IDS[$i]}"
            tr_feito "$id" || continue
            if [ "$tipo" = form ]; then
                form_feito "$id" && continue
                pega_trava "$d/.form.trava" || continue
                diz "FORM $id ${OPCS[$i]} $(instantes "$i")"
                (
                    # shellcheck disable=SC2046,SC2086
                    if python3 -I "$PY/form_kl.py" "$EXE" "$d/form" prog=chuva gradiente=modal M="$M" \
                        iteracoes="$FORM_IT" paralelos=1 threads="$FORM_THREADS" $BASE ${OPCS[$i]} $(instantes "$i") \
                        pressoes="$d/pressoes.txt" nystrom="$NYS_LOCAL" > "$d/form.log" 2>&1; then
                        diz "FORM $id pronto: $(grep -m 1 'beta = ' "$d/form.log" | tr -s ' ')"
                    else
                        echo "ERRO $(date '+%F %T')" >> "$d/.form.trava"
                        diz "FORM $id FALHOU (veja $d/form.log)"
                    fi
                ) &
                n=$((n + FORM_THREADS))
            else
                det_feito "$id" && continue
                pega_trava "$d/.det.trava" || continue
                (
                    mkdir -p "$d/det" && cd "$d/det" || exit 1
                    # shellcheck disable=SC2086
                    "$EXE" n=0 det=1 $BASE ${OPCS[$i]} tempos=$TODOS M="$M" nystrom="$NYS_LOCAL" \
                        pressoes="$d/pressoes.txt" saida=det > det.log 2>&1
                    if det_feito "$id"; then
                        diz "deterministico $id: FS $(awk -F, 'NR > 1 { printf "%s h %s; ", $2, $3 }' det_0_-1.csv)"
                    else
                        echo "ERRO $(date '+%F %T')" >> "$d/.det.trava"
                        diz "deterministico $id FALHOU (veja $d/det/det.log)"
                    fi
                ) &
                n=$((n + 1))
            fi
        done
    done
    wait
}

# fase 2: o line sampling do caso <i> (com a trava ja pega)
ls_caso() {
    local i="$1" id="${IDS[$1]}" d="$SAIDA/${IDS[$1]}" pid rc
    rm -f "$d/ls/PARAR"
    base_local
    diz "line sampling $id"
    python3 -I "$PY/line_sampling.py" "$EXE" "$d/form/form.json" "$d/ls" K="$M" semente=7 cov="$COV" nmin="$NMIN" \
        n="${NLS[$i]}" paralelos="$PARALELO" threads=0 nystrom="$NYS_LOCAL" pressoes="$d/pressoes.txt" >> "$d/ls.log" 2>&1 &
    pid=$!
    while kill -0 "$pid" 2> /dev/null; do
        [ -e "$SAIDA/PARAR" ] && [ ! -e "$d/ls/PARAR" ] && touch "$d/ls/PARAR"
        sleep 30
    done
    wait "$pid"
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "ERRO $(date '+%F %T')" >> "$d/.ls.trava"
        diz "line sampling $id FALHOU (codigo $rc; veja $d/ls.log)"
        return
    fi
    if [ -s "$d/ls/ls.json" ] && ! grep -q '"motivo": "PARAR' "$d/ls/ls.json"; then
        touch "$d/FEITO_LS"
        diz "line sampling $id pronto: $(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1]))
print("P_f %.4f%%, CoV %.1f%%, %d linhas (%d perdidas), %d analises" % (100*d["P_f"], 100*(d["CoV"] or 0), d["linhas"], d["perdidas"], d["analises"]) if d["P_f"] is not None else "sem P_f")' "$d/ls/ls.json")"
    else
        diz "line sampling $id parado (PARAR)"
    fi
    rm -f "$d/.ls.trava" # feito, ou parado (outro job retoma depois)
}

# o primeiro caso com o line sampling livre (FORM pronto, sem FEITO_LS e sem trava)
livre_ls() {
    local i
    for i in "${!IDS[@]}"; do
        ls_feito "${IDS[$i]}" && continue
        form_feito "${IDS[$i]}" || continue
        [ -e "$SAIDA/${IDS[$i]}/.ls.trava" ] && continue
        echo "$i"
        return 0
    done
    return 1
}

# alguma tarefa das fases 0 e 1 ainda por fazer (rodando em outro job ou sem dono)?
pendente() {
    local id
    base_feita || { com_erro "$BASE_NYS.trava" || return 0; }
    for id in "${IDS[@]}"; do
        tr_feito "$id" || { com_erro "$SAIDA/$id/.tr.trava" || return 0; }
        form_feito "$id" || { com_erro "$SAIDA/$id/.form.trava" || return 0; }
    done
    return 1
}

# fase 3: os blocos do Monte Carlo (e da conferencia mc13), intercalados entre os casos: "<i> <pasta> <a> <b>"
mc_tarefas() {
    local i a b k nb maior=0
    for i in "${!IDS[@]}"; do [ "${NMCS[$i]}" -gt "$maior" ] && maior="${NMCS[$i]}"; done
    for i in "${!IDS[@]}"; do # a conferencia primeiro: e pequena e confirma os instantes
        for ((a = 0; a < NCHECK; a += BLOCO_MC)); do
            b=$((a + BLOCO_MC - 1)); [ "$b" -ge "$NCHECK" ] && b=$((NCHECK - 1))
            echo "$i mc13 $a $b"
        done
    done
    nb=$(((maior + BLOCO_MC - 1) / BLOCO_MC))
    for ((k = 0; k < nb; k++)); do
        a=$((k * BLOCO_MC))
        for i in "${!IDS[@]}"; do
            [ "$a" -lt "${NMCS[$i]}" ] || continue
            b=$((a + BLOCO_MC - 1)); [ "$b" -ge "${NMCS[$i]}" ] && b=$((NMCS[$i] - 1))
            echo "$i mc $a $b"
        done
    done
}
# pega ate PARALELO blocos livres (com o transiente e a base prontos), roda todos e espera; 1 se nao havia bloco livre
mc_onda() {
    local n=0 i p a b d ni ins
    base_feita || return 1
    base_local
    while read -r i p a b; do
        [ "$n" -ge "$PARALELO" ] && break
        tr_feito "${IDS[$i]}" || continue
        d="$SAIDA/${IDS[$i]}/$p"
        if [ "$p" = mc13 ]; then ins="tempos=$TODOS" ni=$NTODOS; else ins="$(instantes "$i")" ni=$(($(tr -cd , <<< "${ins%% *}" | wc -c) + 1)); fi
        mc_feito "$d" "$a" "$b" "$ni" && continue
        pega_trava "$d/.tarefas/${a}_${b}" || continue
        rm -f "$d/mc_${a}_${b}.csv" "$d/mc_mapa_${a}_${b}.csv"
        (
            cd "$d" || exit 1
            # shellcheck disable=SC2086
            "$EXE" n=$((b - a + 1)) primeira="$a" $BASE ${OPCS[$i]} $ins M="$M" nystrom="$NYS_LOCAL" \
                pressoes="$SAIDA/${IDS[$i]}/pressoes.txt" saida=mc > "log_${a}_${b}.txt" 2>&1 ||
                echo "ERRO $(date '+%F %T')" >> ".tarefas/${a}_${b}"
        ) &
        n=$((n + 1))
    done < <(mc_tarefas)
    [ "$n" -gt 0 ] || return 1
    wait
    return 0
}

diz "inicio em $(hostname -s): ${#IDS[@]} casos, k_s $KS, M $M, transiente $TRANS, instantes $TEMPOS, PARALELO" \
    "$PARALELO, CoV alvo $COV, linhas $NMIN a $N, FORM ate $FORM_IT iteracoes, MC ate $NMC_MAX por caso (blocos de" \
    "$BLOCO_MC), conferencia $NCHECK, executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
while :; do
    [ -e "$SAIDA/PARAR" ] && { diz "achei $SAIDA/PARAR: saindo"; break; }
    fase0
    fase1
    if i="$(livre_ls)"; then
        pega_trava "$SAIDA/${IDS[$i]}/.ls.trava" || continue
        ls_caso "$i"
        continue
    fi
    if pendente; then # nenhum line sampling livre ainda: adianta o Monte Carlo (ou espera um pouco)
        mc_onda || sleep 60
        continue
    fi
    mc_onda && continue # os line samplings estao feitos ou com outros jobs: o Monte Carlo
    diz "nada mais livre: saindo"
    break
done
diz "fim"
