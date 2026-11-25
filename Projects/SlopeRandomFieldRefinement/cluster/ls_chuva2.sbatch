#!/bin/bash
# ls_chuva2.sbatch: FORM modal + line sampling da campanha chuva2 (a "Secao 6 da chuva": talude de Rahardjo, variante iii,
# esqueleto rigido e ar passivo, SlopeRainMonteCarlo), um LS por linha de casos_chuva2_ls.txt, com as opcoes do caso do
# Monte Carlo (casos_chuva2.txt, as mesmas do job.sbatch CAMPANHA=chuva2) e a falha nos instantes T da linha:
#   g = ln(FS_min), FS_min = o menor FS dos instantes de T (tempos=T e adapta=T no lugar dos do caso); falha: FS_min < 1
#   (primeira passagem nesses instantes; um instante sem equilibrio no SRM e falha, FS = 0,005; um transiente que parou
#   tira a linha da conta). Referencia: T = 24 h (o fim da chuva; o menor FS do f50 cai em 24 h em 99,3% das realizacoes).
# Direcao: o FORM com o gradiente modal em c' e phi' (form_kl.py prog=chuva gradiente=modal, M termos por campo, 1
# analise por iteracao, FORM_THREADS threads): o gradiente modal da 0 em k (o fator congela as pressoes), entao a direcao
# nao tem componente em k e o k_s e aleatorio ENTRE as linhas (na parte ortogonal), nao ao longo delas. O LS nao tem vies
# para nenhuma direcao. Line sampling: line_sampling.py busca=cerca (raiz cercada pelo sinal de FS_min - 1, refino ITP),
# K = M termos (o campo do MC: M = 2000 da KL analitica, ou o M= do caso), semente 7, PARALELO linhas em voo, ate
# CoV(P_f) <= COV (no minimo NMIN e no maximo N_ls linhas, limitado a N).
# REUSO=1 (padrao 0): o transiente calculado uma vez por linha (e uma vez no FORM): ao longo da linha so c' e phi' mudam,
# entao a 1a analise grava as pressoes (grava_pressoes=) e as outras so calculam o fator (pressoes=); o FS e o mesmo bit
# a bit. Precisa do executavel com o patch do campo de k (a chave do transiente com o campo de k; o submeter confere).
# Os jobs (submeter_ls_chuva2.sh) dividem o trabalho por travas em $SAIDA (noclobber, como no hu2024_casos.sbatch):
#   1. o FORM e o deterministico (as medias, det=1, com os instantes do caso) de cada LS: ate PARALELO nucleos por job
#      (FORM_THREADS por FORM, 1 por deterministico), todos ao mesmo tempo, e espera;
#   2. o line sampling de cada LS com o FORM pronto, um por vez por job (na ordem do arquivo: a referencia primeiro).
#   Enquanto ha FORM rodando em outro job e nenhum LS livre, o job espera; sem nada livre, sai.
# Retomada: o FORM interrompido recomeca do zero; o line sampling continua do linhas.csv. As travas de um job que morreu
# sao liberadas pelo submeter_ls_chuva2.sh quando nao ha jobs lsChuva2 na fila. Parar sem perder nada: touch $SAIDA/PARAR
# (o LS para de submeter linhas e termina as em voo); para voltar: rm $SAIDA/PARAR e ./submeter_ls_chuva2.sh de novo.
# Saidas por LS em $SAIDA/<id>/: det/ (det_0_-1.csv), form/ (form.json), form.log, ls/ (ls.json, linhas.csv, analises.csv,
# perdidas.csv; custos.csv com REUSO=1), ls.log, FEITO_LS e as travas. Eventos em $SAIDA/eventos.log. Resumo contra o MC
# direto da chuva2 (o mesmo evento): resumo_ls_chuva2.py.
# Uso: ./submeter_ls_chuva2.sh [nos]. Ambiente: SAIDA ($DATA/chuva2_ls), SAIDA_MC ($DATA/chuva2: o @SAIDA@ das opcoes dos
# casos), IREF (31.05), PARALELO (24), COV (0.05), NMIN (30), N (600), TOL (0.01, em c), FORM_IT (15), FORM_THREADS (4),
# FORM_EXTRA (opcoes a mais do form_kl.py, ex. algoritmo=ihlrf), REUSO (0), LS_CASOS (todos, ou ids separados por
# virgulas), CASOS_TXT, LS_TXT, EXE (bin_chuva2/SlopeRainMonteCarlo; o do MC: EXE=$PWD/bin/SlopeRainMonteCarlo).
#SBATCH --job-name=lsChuva2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}" || exit 1
[ -n "${SLURM_JOB_ID:-}" ] && source ../../LearningMohrCoulomb2/cluster/modulos.sh
DADOS="${DATA:-/mnt/data/${USER:-$(id -un)}}"
SAIDA="${SAIDA:-$DADOS/chuva2_ls}" SAIDA_MC="${SAIDA_MC:-$DADOS/chuva2}" IREF="${IREF:-31.05}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}"
COV="${COV:-0.05}" NMIN="${NMIN:-30}" N="${N:-600}" TOL="${TOL:-0.01}"
FORM_IT="${FORM_IT:-15}" FORM_THREADS="${FORM_THREADS:-4}" FORM_EXTRA="${FORM_EXTRA:-}" REUSO="${REUSO:-0}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_chuva2.txt}" LS_TXT="${LS_TXT:-$PWD/casos_chuva2_ls.txt}" LS_CASOS="${LS_CASOS:-todos}"
EXE="${EXE:-$PWD/bin_chuva2/SlopeRainMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)"
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_ls_chuva2.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }
OPC_REUSO=""
[ "$REUSO" = 1 ] && OPC_REUSO="reuso=1"

# os casos do MC: id -> opcoes, como o tarefas.sh as passa (IREF e @SAIDA@ trocados, sem as opcoes da fila)
declare -A OPC_MC
while read -r id fase n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    opc="${opc//IREF/$IREF}" resto=""
    set -f
    for w in $opc; do
        case "$w" in prog=* | espera=* | ref=* | par=* | cov_alvo=* | cov_n=* | cov_falhas=*) ;; *) resto+="${resto:+ }$w" ;; esac
    done
    set +f
    OPC_MC[$id]="${resto//@SAIDA@/$SAIDA_MC}"
done < "$CASOS_TXT"
# os LS: id, opcoes (as do caso com tempos=T), T, teto de linhas e M (LS_CASOS: so os ids pedidos)
LIDS=() LCASO=() LOPC=() LT=() LN=() LM=()
while read -r id caso t nls; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$LS_CASOS" != todos ] && [[ ",$LS_CASOS," != *",$id,"* ]]; then continue; fi
    [ -n "${OPC_MC[$caso]+x}" ] || { echo "LS $id: o caso $caso nao esta em $CASOS_TXT" >&2; exit 1; }
    opc="" m=2000
    set -f
    for w in ${OPC_MC[$caso]}; do
        case "$w" in
        tempos=* | adapta=*) ;;
        piloto=* | desloc_xi=* | banda=*) echo "LS $id: o caso $caso e de amostragem por importancia ($w)" >&2; exit 1 ;;
        M=*) m="${w#M=}" opc+="${opc:+ }$w" ;;
        *) opc+="${opc:+ }$w" ;;
        esac
    done
    set +f
    [ "$nls" -gt "$N" ] && nls="$N"
    LIDS+=("$id") LCASO+=("$caso") LOPC+=("$opc tempos=$t") LT+=("$t") LN+=("$nls") LM+=("$m")
done < "$LS_TXT"
[ ${#LIDS[@]} -gt 0 ] || { echo "nenhum LS (LS_CASOS=$LS_CASOS em $LS_TXT)" >&2; exit 1; }

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
det_feito() { [ -s "$SAIDA/$1/det/det_0_-1.csv" ] && [ "$(wc -l < "$SAIDA/$1/det/det_0_-1.csv")" -gt 1 ]; }
ls_feito() { [ -e "$SAIDA/$1/FEITO_LS" ]; }

# a base de Nystrom (se o caso tem nystrom=) no no: copia local, ou a original se a copia falhar; as opcoes <i> com ela
nys_local() {
    local b="$1" loc
    loc="${TMPDIR:-/tmp}/chuva2ls_$(basename "$b")"
    if [ ! -s "$loc" ] || [ "$(stat -c %s "$loc")" != "$(stat -c %s "$b")" ]; then
        cp "$b" "$loc.$JOB" 2> /dev/null && mv -f "$loc.$JOB" "$loc" || loc="$b"
    fi
    echo "$loc"
}
opcoes() { # <opcoes>: com a base de Nystrom local no lugar da original
    local w out=""
    set -f
    for w in $1; do
        case "$w" in nystrom=*) w="nystrom=$(nys_local "${w#nystrom=}")" ;; esac
        out+="${out:+ }$w"
    done
    set +f
    echo "$out"
}

# fase 1: FORM e deterministico (com os instantes todos do caso); pega ate PARALELO nucleos de tarefas livres e espera
fase1() {
    local tipo i id d n=0 opc
    for tipo in form det; do
        for i in "${!LIDS[@]}"; do
            [ "$n" -ge "$PARALELO" ] && break 2
            id="${LIDS[$i]}" d="$SAIDA/${LIDS[$i]}"
            if [ "$tipo" = form ]; then
                form_feito "$id" && continue
                pega_trava "$d/.form.trava" || continue
                opc="$(opcoes "${LOPC[$i]}")"
                diz "FORM $id: caso ${LCASO[$i]}, $opc${OPC_REUSO:+ (reuso)}"
                (
                    # shellcheck disable=SC2086
                    if python3 -I "$PY/form_kl.py" "$EXE" "$d/form" prog=chuva gradiente=modal M="${LM[$i]}" \
                        iteracoes="$FORM_IT" paralelos=1 threads="$FORM_THREADS" $FORM_EXTRA $OPC_REUSO $opc \
                        > "$d/form.log" 2>&1; then
                        diz "FORM $id pronto: $(grep -m 1 'beta = ' "$d/form.log" | tr -s ' ')$(grep -q 'NAO convergiu' "$d/form.log" && echo ' (nao convergiu)')"
                    else
                        echo "ERRO $(date '+%F %T')" >> "$d/.form.trava"
                        diz "FORM $id FALHOU (veja $d/form.log)"
                    fi
                ) &
                n=$((n + FORM_THREADS))
            else
                det_feito "$id" && continue
                pega_trava "$d/.det.trava" || continue
                opc="$(opcoes "${OPC_MC[${LCASO[$i]}]}")"
                (
                    mkdir -p "$d/det" && cd "$d/det" || exit 1
                    # shellcheck disable=SC2086
                    "$EXE" n=0 det=1 $opc saida=det > det.log 2>&1
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

# fase 2: o line sampling do LS <i> (com a trava ja pega)
ls_item() {
    local i="$1" id="${LIDS[$1]}" d="$SAIDA/${LIDS[$1]}" pid rc nys
    rm -f "$d/ls/PARAR"
    nys="$(opcoes "${LOPC[$i]}" | tr ' ' '\n' | grep '^nystrom=')"
    diz "line sampling $id (ate ${LN[$i]} linhas, CoV alvo $COV, K ${LM[$i]}${OPC_REUSO:+, reuso})"
    # shellcheck disable=SC2086
    python3 -I "$PY/line_sampling.py" "$EXE" "$d/form/form.json" "$d/ls" K="${LM[$i]}" semente=7 cov="$COV" \
        nmin="$NMIN" n="${LN[$i]}" paralelos="$PARALELO" threads=0 busca=cerca tol="$TOL" $OPC_REUSO $nys \
        >> "$d/ls.log" 2>&1 &
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

# o primeiro LS livre (FORM pronto, sem FEITO_LS e sem trava), na ordem do arquivo
livre_ls() {
    local i
    for i in "${!LIDS[@]}"; do
        ls_feito "${LIDS[$i]}" && continue
        form_feito "${LIDS[$i]}" || continue
        [ -e "$SAIDA/${LIDS[$i]}/.ls.trava" ] && continue
        echo "$i"
        return 0
    done
    return 1
}
# algum FORM ainda por fazer (rodando em outro job), sem ERRO?
pendente() {
    local id
    for id in "${LIDS[@]}"; do
        form_feito "$id" || { com_erro "$SAIDA/$id/.form.trava" || return 0; }
    done
    return 1
}

diz "inicio em $(hostname -s): ${#LIDS[@]} LS (${LIDS[*]}), PARALELO $PARALELO, CoV alvo $COV, linhas $NMIN a $N, tol" \
    "$TOL, FORM ate $FORM_IT iteracoes com $FORM_THREADS threads${FORM_EXTRA:+ ($FORM_EXTRA)}, reuso $REUSO, IREF $IREF," \
    "executavel $EXE ($(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null))"
while :; do
    [ -e "$SAIDA/PARAR" ] && { diz "achei $SAIDA/PARAR: saindo"; break; }
    fase1
    if i="$(livre_ls)"; then
        pega_trava "$SAIDA/${LIDS[$i]}/.ls.trava" || continue
        ls_item "$i"
        continue
    fi
    if pendente; then # um FORM ainda roda em outro job: espera para fazer o LS dele
        sleep 60
        continue
    fi
    diz "nada mais livre: saindo"
    break
done
diz "fim"
