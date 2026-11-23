#!/bin/bash
# liu2024_casos.sbatch: o caso B, Liu, Han, Huang, Huang e Cheng (2024) (talude de Griffiths e Fenton, 5 m a 1:2, com a
# percolacao permanente de q = 1e-7 m/s; k_s, c' e phi' lognormais com o nucleo gaussiano 30/3 m, r(c, phi) = -0,5 e
# phi^b = 2 phi'/3 por ponto), com o SlopeRainMonteCarlo: FORM com o gradiente modal, line sampling e Monte Carlo direto
# por cenario (casos_liu2024.txt), como o hu2024_casos.sbatch. Num cenario com k_s fixo (covk=0 nas opcoes: o cenario 2)
# a solucao permanente e a mesma em todas as realizacoes: ela e calculada uma vez (so_transiente=1 grava_pressoes=) e
# todas as analises a leem (pressoes=); com k_s aleatorio (cenario 3) cada analise faz a sua.
# Os jobs (submeter_liu2024.sh) dividem o trabalho por travas em $SAIDA (noclobber, como no tarefas.sh), em fases:
#   0. a base de Nystrom do talude (M termos, nucleo gauss 30/3 m, grade 0,2 x 0,1 m; nystrom_base.py com o python do
#      venv) e a solucao permanente dos cenarios com k_s fixo: ate PARALELO tarefas por job;
#   1. o deterministico (as medias) e o FORM (M termos por campo, gradiente modal, ate FORM_IT iteracoes de 1 analise
#      com FORM_THREADS threads) de cada cenario; os termos de k saem 0 do gradiente modal (o fator congela as pressoes):
#      a direcao do FORM fica so em c' e phi', e o LS continua sem vies;
#   2. o line sampling de cada cenario com o FORM pronto, um por vez com PARALELO linhas em voo, ate CoV(P_f) <= COV (no
#      minimo NMIN e no maximo o teto do cenario, limitado a N linhas), semente 7;
#   3. o Monte Carlo direto: min(N_mc do cenario, NMC_MAX) realizacoes por cenario em blocos de BLOCO_MC, intercalados
#      entre os cenarios (a realizacao n tem os mesmos numeros aleatorios em todos).
# O FS e o do SRM em t = 0 (o permanente, sem chuva depois: I=0, tempos=0); a falha e FS < 1.
# Retomada, parada (touch $SAIDA/PARAR) e liberacao das travas: como no hu2024_casos.sbatch.
# Saidas por cenario em $SAIDA/<id>/: pressoes.txt (so com k_s fixo), det/, form/, form.log, ls/, ls.log, FEITO_LS, mc/ e as
# travas. Eventos em $SAIDA/eventos.log. Resumo contra o artigo: resumo_liu2024.py.
# Uso: ./submeter_liu2024.sh [nos]. Ambiente: SAIDA ($DATA/liu2024), M (100), BASE_NYS ($DATA/liu2024/base_liu_M$M.txt),
# TRANS (tri3), NREF (3: os niveis do fator; o 4o muda o FS em ~0,3% e dobra o custo), PARALELO (24), COV (0.05), NMIN
# (30), N (600), FORM_IT (10), FORM_THREADS (4), NMC_MAX (10000), BLOCO_MC (10), CASOS (todos, ou ids separados por
# virgulas), EXE (bin_liu/SlopeRainMonteCarlo), PYV ($HOME/venv_la/bin/python).
#SBATCH --job-name=liuB
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
M="${M:-100}"
SAIDA="${SAIDA:-$DADOS/liu2024}"
BASE_NYS="${BASE_NYS:-$DADOS/liu2024/base_liu_M$M.txt}"
TRANS="${TRANS:-tri3}" NREF="${NREF:-3}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}" FORM_IT="${FORM_IT:-10}" FORM_THREADS="${FORM_THREADS:-4}"
NMC_MAX="${NMC_MAX:-10000}" BLOCO_MC="${BLOCO_MC:-10}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_liu2024.txt}" CASOS="${CASOS:-todos}"
EXE="${EXE:-$PWD/bin_liu/SlopeRainMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)" PYV="${PYV:-$HOME/venv_la/bin/python}"
TALUDE=0,0,30,0,30,5,20,5,10,10,0,10
# o caso B: geometria, nivel inicial reto de (0; 6) a (30; 1) (o das cargas laterais de 6 e 1 m; so o chute do
# permanente), Fredlund-Xing a = 5 kPa, n = 2, m = 1, theta_s = 0,4, k pela integral de Fredlund (psi_aev = 2,6 kPa, a
# tangente; psi_r = 3000 kPa), k_s = 1e-5 m/s, q = 1e-7 m/s (0,36 mm/h) por metro inclinado (o fluxo unitario do
# GeoStudio), chi = tan(2 phi'/3)/tan(phi') por ponto, gamma = 20, c' = 10 kPa (CoV 0,3), phi' = 15 graus (CoV 0,2),
# r(c, phi) = -0,5 (Tabela 1); as malhas. O CoV de k_s (0,6 ou 0) vai nas opcoes do cenario
BASE="geometria=$TALUDE na=2.666666667 na_beta=9.462322208 carga_esq=6 carga_dir=1 fx=5,2,1,0.4 kfredlund=2.6 ks=1e-5"
BASE+=" q0=0.36 I=0 chuva=0 tempos=0 chi=ponto phib_razao=0.666666667 c=10 phi=15 gamma=20 covc=0.3 covphi=0.2"
BASE+=" rcphi=-0.5 transiente=$TRANS fator=tri1 nref=$NREF"
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_hu2024.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }

# os cenarios: id, P_f do artigo, N do MC, teto de linhas do LS e opcoes (CASOS: so os ids pedidos); REUSA: 1 se
# covk=0 (a solucao permanente e calculada uma vez e lida por todas as analises)
IDS=() OPCS=() NMCS=() NLS=() REUSA=()
while read -r id pf nmc nls opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ] && [[ ",$CASOS," != *",$id,"* ]]; then continue; fi
    [ "$nmc" -gt "$NMC_MAX" ] && nmc="$NMC_MAX"
    [ "$nls" -gt "$N" ] && nls="$N"
    r=0; [[ " $opc " == *" covk=0 "* ]] && r=1
    IDS+=("$id") OPCS+=("$opc") NMCS+=("$nmc") NLS+=("$nls") REUSA+=("$r")
done < "$CASOS_TXT"
[ ${#IDS[@]} -gt 0 ] || { echo "nenhum caso (CASOS=$CASOS em $CASOS_TXT)" >&2; exit 1; }
# pressoes= do cenario <i> (so com k_s fixo)
press() { [ "${REUSA[$1]}" = 1 ] && echo "pressoes=$SAIDA/${IDS[$1]}/pressoes.txt"; }
indice() { local i; for i in "${!IDS[@]}"; do [ "${IDS[$i]}" = "$1" ] && echo "$i" && return; done; }

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
tr_feito() { [ "${REUSA[$(indice "$1")]}" = 0 ] || [ -s "$SAIDA/$1/pressoes.txt" ]; }
form_feito() { [ -s "$SAIDA/$1/form/form.json" ]; }
det_feito() { [ -s "$SAIDA/$1/det/det_0_-1.csv" ] && [ "$(wc -l < "$SAIDA/$1/det/det_0_-1.csv")" -gt 1 ]; }
ls_feito() { [ -e "$SAIDA/$1/FEITO_LS" ]; }
mc_feito() { # <pasta mc> <a> <b> <instantes por realizacao>: o csv do bloco com as (b - a + 1) x instantes linhas
    local f="$1/mc_$2_$3.csv"
    [ -s "$f" ] && [ $(($(wc -l < "$f") - 1)) -ge $((($3 - $2 + 1) * $4)) ]
}

# a base de Nystrom no no: copia local (todas as analises do no leem a mesma), ou a de $DATA se a copia falhar
NYS_LOCAL=""
base_local() {
    local loc="${TMPDIR:-/tmp}/liu2024_$(basename "$BASE_NYS")"
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
            if "$PYV" -I "$PY/nystrom_base.py" "$BASE_NYS.tmp" M="$M" talude="$TALUDE" nucleo=gauss lh=30 lv=3 dx=0.2 dy=0.1 \
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
        diz "permanente $id (${OPCS[$i]})"
        (
            # shellcheck disable=SC2086
            if PZ_THREADS=4 "$EXE" $BASE ${OPCS[$i]} so_transiente=1 grava_pressoes="$d/pressoes.txt" \
                saida="$d/tr" > "$d/transiente.log" 2>&1; then
                diz "permanente $id pronto: $(grep -m 1 'passos' "$d/transiente.log" | tr -s ' ')"
            else
                echo "ERRO $(date '+%F %T')" >> "$d/.tr.trava"
                diz "permanente $id FALHOU (veja $d/transiente.log)"
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
                diz "FORM $id ${OPCS[$i]}"
                (
                    # shellcheck disable=SC2046,SC2086
                    if python3 -I "$PY/form_kl.py" "$EXE" "$d/form" prog=chuva gradiente=modal M="$M" \
                        iteracoes="$FORM_IT" paralelos=1 threads="$FORM_THREADS" $BASE ${OPCS[$i]} $(press "$i") \
                        nystrom="$NYS_LOCAL" > "$d/form.log" 2>&1; then
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
                    "$EXE" n=0 det=1 $BASE ${OPCS[$i]} M="$M" nystrom="$NYS_LOCAL" $(press "$i") saida=det > det.log 2>&1
                    if det_feito "$id"; then
                        diz "deterministico $id: FS $(awk -F, 'NR > 1 { printf "%s", $3 }' det_0_-1.csv)"
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
        n="${NLS[$i]}" paralelos="$PARALELO" threads=0 nystrom="$NYS_LOCAL" $(press "$i") >> "$d/ls.log" 2>&1 &
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

# fase 3: os blocos do Monte Carlo, intercalados entre os cenarios: "<i> <a> <b>"
mc_tarefas() {
    local i a b k nb maior=0
    for i in "${!IDS[@]}"; do [ "${NMCS[$i]}" -gt "$maior" ] && maior="${NMCS[$i]}"; done
    nb=$(((maior + BLOCO_MC - 1) / BLOCO_MC))
    for ((k = 0; k < nb; k++)); do
        a=$((k * BLOCO_MC))
        for i in "${!IDS[@]}"; do
            [ "$a" -lt "${NMCS[$i]}" ] || continue
            b=$((a + BLOCO_MC - 1)); [ "$b" -ge "${NMCS[$i]}" ] && b=$((NMCS[$i] - 1))
            echo "$i $a $b"
        done
    done
}
# pega ate PARALELO blocos livres (com o permanente e a base prontos), roda todos e espera; 1 se nao havia bloco livre
mc_onda() {
    local n=0 i a b d
    base_feita || return 1
    base_local
    while read -r i a b; do
        [ "$n" -ge "$PARALELO" ] && break
        tr_feito "${IDS[$i]}" || continue
        d="$SAIDA/${IDS[$i]}/mc"
        mc_feito "$d" "$a" "$b" 1 && continue
        pega_trava "$d/.tarefas/${a}_${b}" || continue
        rm -f "$d/mc_${a}_${b}.csv" "$d/mc_mapa_${a}_${b}.csv"
        (
            cd "$d" || exit 1
            # shellcheck disable=SC2086
            "$EXE" n=$((b - a + 1)) primeira="$a" $BASE ${OPCS[$i]} M="$M" nystrom="$NYS_LOCAL" $(press "$i") saida=mc \
                > "log_${a}_${b}.txt" 2>&1 || echo "ERRO $(date '+%F %T')" >> ".tarefas/${a}_${b}"
        ) &
        n=$((n + 1))
    done < <(mc_tarefas)
    [ "$n" -gt 0 ] || return 1
    wait
    return 0
}

diz "inicio em $(hostname -s): ${#IDS[@]} cenarios, M $M, malha do permanente $TRANS, PARALELO $PARALELO, CoV alvo" \
    "$COV, linhas $NMIN a $N, FORM ate $FORM_IT iteracoes, MC ate $NMC_MAX por cenario (blocos de $BLOCO_MC)," \
    "executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
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
