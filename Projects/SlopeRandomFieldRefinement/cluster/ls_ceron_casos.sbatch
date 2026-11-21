#!/bin/bash
# ls_ceron_casos.sbatch: as Secoes 6.1 a 6.5 de Vargas Ceron et al. (2025) com a agua do artigo (contorno=ceron agua=sat
# trebaix=0.0001: percolacao saturada depois do rebaixamento instantaneo, Eqs. 19-21), por FORM com o gradiente modal e
# line sampling (como o ls_ceron.sbatch, que roda so a referencia), mais um Monte Carlo direto pareado para as
# densidades, a media e o CoV de Gamma das figuras do artigo. Os 26 casos estao em casos_ls_ceron.txt; os jobs
# (submeter_ls_ceron.sh) dividem o trabalho por travas em $SAIDA (noclobber, como no tarefas.sh), em tres fases:
#   1. o FORM (M = 2000 termos por campo, ate FORM_IT iteracoes de 1 analise) e o deterministico (as medias, det=1) de
#      cada caso: cada job pega ate PARALELO dessas tarefas (os FORM primeiro), roda todas ao mesmo tempo e espera;
#   2. o line sampling de cada caso com o FORM pronto, um caso por vez com PARALELO linhas em voo, ate CoV(P_f) <= COV
#      (no minimo NMIN e no maximo N linhas), semente 7 em todos os casos (os mesmos Z: numeros aleatorios comuns);
#   3. o Monte Carlo direto: NMC realizacoes por caso (padrao 500), em blocos de BLOCO_MC, intercalados entre os casos
#      (o bloco 0 de todos, depois o 1, ...); a realizacao n tem os mesmos numeros aleatorios em todos os casos (a
#      semente padrao do SlopeDrawdownMonteCarlo), como no artigo. Enquanto ha FORM rodando em outro job e nenhum line
#      sampling livre, o job adianta o Monte Carlo em vez de esperar.
# Base comum dos casos: contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5 (Gamma no permanente) transiente=tri3
# fator=tri1 nref=5 (com 3 niveis o Gamma fica ~4% alto); o caso acrescenta as suas opcoes (covc=, s=, alfa=, hw=...).
# Retomada: o FORM interrompido recomeca do zero; o line sampling continua do linhas.csv; o bloco do MC incompleto e
# refeito. As travas de um job que morreu sao liberadas pelo submeter_ls_ceron.sh quando nao ha jobs lsCeronC na fila.
# Parar sem perder nada: touch $SAIDA/PARAR (os line samplings em curso param de submeter linhas e terminam as em voo;
# os blocos do MC em curso terminam). Para voltar: rm $SAIDA/PARAR e ./submeter_ls_ceron.sh de novo.
# Saidas por caso em $SAIDA/<id>/: det/ (det_0_-1.csv), form/ (form.json), form.log, ls/ (ls.json, linhas.csv,
# analises.csv, perdidas.csv), ls.log, FEITO_LS, mc/ (mc_<a>_<b>.csv) e as travas (.form.trava, .det.trava,
# .ls.trava, mc/.tarefas/). Eventos de todos os jobs em $SAIDA/eventos.log. Resumo: resumo_ls_ceron.py.
# Uso: ./submeter_ls_ceron.sh [nos]  (ou, para um job so: sbatch ls_ceron_casos.sbatch). Ambiente: SAIDA
# ($DATA/ceron_ls_casos), PARALELO (24), COV (0.05), NMIN (30), N (600), FORM_IT (15), NMC (500), BLOCO_MC (5),
# CASOS (todos, ou ids separados por virgulas), EXE (bin_ceron/SlopeDrawdownMonteCarlo).
#SBATCH --job-name=lsCeronC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/ceron_ls_casos}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}" FORM_IT="${FORM_IT:-15}" NMC="${NMC:-500}" BLOCO_MC="${BLOCO_MC:-5}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_ls_ceron.txt}" CASOS="${CASOS:-todos}"
EXE="${EXE:-$PWD/bin_ceron/SlopeDrawdownMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)"
BASE="contorno=ceron agua=sat trebaix=0.0001 tfim=20 dt=5 transiente=tri3 fator=tri1 nref=5"
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_ceron.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }

# os casos: id e opcoes (CASOS: so os ids pedidos)
IDS=() OPCS=()
while read -r id secao pf cvpf gdet opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ] && [[ ",$CASOS," != *",$id,"* ]]; then continue; fi
    [ "$opc" = "-" ] && opc=""
    IDS+=("$id") OPCS+=("$opc")
done < "$CASOS_TXT"
[ ${#IDS[@]} -gt 0 ] || { echo "nenhum caso (CASOS=$CASOS em $CASOS_TXT)" >&2; exit 1; }

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
form_feito() { [ -s "$SAIDA/$1/form/form.json" ]; }
det_feito() { [ -s "$SAIDA/$1/det/det_0_-1.csv" ]; }
ls_feito() { [ -e "$SAIDA/$1/FEITO_LS" ]; }
mc_feito() { # <indice do caso> <a> <b>: o csv do bloco com as b - a + 1 realizacoes
    local f="$SAIDA/${IDS[$1]}/mc/mc_$2_$3.csv"
    [ -s "$f" ] && [ $(($(wc -l < "$f") - 1)) -ge $(($3 - $2 + 1)) ]
}

# fase 1: FORM e deterministico; pega ate PARALELO tarefas livres, roda todas e espera
fase1() {
    local tipo i id d n=0
    for tipo in form det; do
        for i in "${!IDS[@]}"; do
            [ "$n" -ge "$PARALELO" ] && break 2
            id="${IDS[$i]}" d="$SAIDA/${IDS[$i]}"
            if [ "$tipo" = form ]; then
                form_feito "$id" && continue
                pega_trava "$d/.form.trava" || continue
                diz "FORM $id ${OPCS[$i]}"
                (
                    # shellcheck disable=SC2086
                    if python3 -I "$PY/form_kl.py" "$EXE" "$d/form" prog=rebaixamento gradiente=modal M=2000 \
                        iteracoes="$FORM_IT" paralelos=1 $BASE ${OPCS[$i]} > "$d/form.log" 2>&1; then
                        diz "FORM $id pronto: $(grep -m 1 'beta = ' "$d/form.log" | tr -s ' ')"
                    else
                        echo "ERRO $(date '+%F %T')" >> "$d/.form.trava"
                        diz "FORM $id FALHOU (veja $d/form.log)"
                    fi
                ) &
            else
                det_feito "$id" && continue
                pega_trava "$d/.det.trava" || continue
                (
                    mkdir -p "$d/det" && cd "$d/det" || exit 1
                    # shellcheck disable=SC2086
                    "$EXE" n=0 det=1 $BASE ${OPCS[$i]} saida=det > det.log 2>&1
                    if [ -s det_0_-1.csv ]; then
                        diz "deterministico $id: Gamma $(sed -n 2p det_0_-1.csv | cut -d, -f3,4)"
                    else
                        echo "ERRO $(date '+%F %T')" >> "$d/.det.trava"
                        diz "deterministico $id FALHOU (veja $d/det/det.log)"
                    fi
                ) &
            fi
            n=$((n + 1))
        done
    done
    wait
}

# fase 2: o line sampling do caso <i> (com a trava ja pega)
ls_caso() {
    local i="$1" id="${IDS[$1]}" d="$SAIDA/${IDS[$1]}" pid rc
    rm -f "$d/ls/PARAR"
    diz "line sampling $id"
    python3 -I "$PY/line_sampling.py" "$EXE" "$d/form/form.json" "$d/ls" K=2000 semente=7 cov="$COV" nmin="$NMIN" \
        n="$N" paralelos="$PARALELO" >> "$d/ls.log" 2>&1 &
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
print("P_f %.3f%%, CoV %.1f%%, %d linhas (%d perdidas), %d analises" % (100*d["P_f"], 100*(d["CoV"] or 0), d["linhas"], d["perdidas"], d["analises"]) if d["P_f"] is not None else "sem P_f")' "$d/ls/ls.json")"
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

# algum FORM ainda por fazer (rodando em outro job ou sem dono)?
form_pendente() {
    local id
    for id in "${IDS[@]}"; do
        form_feito "$id" && continue
        com_erro "$SAIDA/$id/.form.trava" && continue
        return 0
    done
    return 1
}

# fase 3: os blocos do Monte Carlo, intercalados entre os casos
mc_tarefas() {
    local nb=$(((NMC + BLOCO_MC - 1) / BLOCO_MC)) k i a b
    for ((k = 0; k < nb; k++)); do
        a=$((k * BLOCO_MC)) b=$((k * BLOCO_MC + BLOCO_MC - 1))
        [ "$b" -ge "$NMC" ] && b=$((NMC - 1))
        for i in "${!IDS[@]}"; do echo "$i $a $b"; done
    done
}
# pega ate PARALELO blocos livres, roda todos e espera; 1 se nao havia bloco livre
mc_onda() {
    local n=0 i a b d
    while read -r i a b; do
        [ "$n" -ge "$PARALELO" ] && break
        mc_feito "$i" "$a" "$b" && continue
        d="$SAIDA/${IDS[$i]}/mc"
        pega_trava "$d/.tarefas/${a}_${b}" || continue
        rm -f "$d/mc_${a}_${b}.csv" "$d/mc_mapa_${a}_${b}.csv"
        (
            cd "$d" || exit 1
            # shellcheck disable=SC2086
            "$EXE" n=$((b - a + 1)) primeira="$a" $BASE ${OPCS[$i]} saida=mc > "log_${a}_${b}.txt" 2>&1 ||
                echo "ERRO $(date '+%F %T')" >> ".tarefas/${a}_${b}"
        ) &
        n=$((n + 1))
    done < <(mc_tarefas)
    [ "$n" -gt 0 ] || return 1
    wait
    return 0
}

diz "inicio em $(hostname -s): ${#IDS[@]} casos, PARALELO $PARALELO, CoV alvo $COV, linhas $NMIN a $N, FORM ate" \
    "$FORM_IT iteracoes, MC $NMC por caso (blocos de $BLOCO_MC), executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
while :; do
    [ -e "$SAIDA/PARAR" ] && { diz "achei $SAIDA/PARAR: saindo"; break; }
    fase1
    if i="$(livre_ls)"; then
        pega_trava "$SAIDA/${IDS[$i]}/.ls.trava" || continue
        ls_caso "$i"
        continue
    fi
    if form_pendente; then # nenhum line sampling livre ainda: adianta o Monte Carlo (ou espera um pouco)
        mc_onda || sleep 60
        continue
    fi
    mc_onda && continue # os line samplings estao feitos ou com outros jobs: o Monte Carlo
    diz "nada mais livre: saindo"
    break
done
diz "fim"
