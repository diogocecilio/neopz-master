#!/bin/bash
# liuwang2023_casos.sbatch: o caso C, Liu e Wang (2023) (talude de 10 m, 74 x 30 m, cargas de 20 e 10 m, van
# Genuchten-Mualem, phi^b = 15 graus fixo, chi = tan(15)/tan(phi'(x)) por ponto; k_s, c' e phi' lognormais com o nucleo
# exponencial 20/2 m por Nystrom; estado inicial permanente com q0 = chuva anual / 8760 h, a chuva anual normal de 1200
# mm e CoV 0,1; chuva (I, D) Pareto generalizada com copula de Frank), com o SlopeRainMonteCarlo, por caso de
# casos_liuwang2023.txt:
#   crpc: o metodo semi-analitico do artigo (Secao 3): a curva critica (CRPC) de cada realizacao (opcao crpc=: 9
#         intensidades, de 100 a 1 mm/h, chuva constante ate 200 h, D_C pela bissecao do FS; as I menores nao sao feitas
#         depois de uma que nao derruba) e P_1 pela integral da densidade de (I, D) acima da curva (crpc_liuwang.py); uma
#         tarefa por realizacao. O caso deterministico (o solo das medias, N = 0) roda uma tarefa por intensidade;
#   mc:   o Monte Carlo direto (chuva_aleatoria=1: I, D e q0 sorteados; FS em 0 e no fim da chuva), em blocos de
#         BLOCO_MC realizacoes.
# As amostras com FS0 < 1 saem da conta (o artigo: 311 de 5484); a P_F1 e por evento de chuva.
# Os jobs (submeter_liuwang2023.sh) dividem o trabalho por travas em $SAIDA (noclobber, como no hu2024_casos.sbatch):
# primeiro a base de Nystrom (uma tarefa), depois ondas de ate PARALELO tarefas (as do caso deterministico e as CRPC
# primeiro, intercaladas entre os casos; depois os blocos do MC). Retomada: a tarefa incompleta e refeita; as travas de
# jobs que morreram sao liberadas pelo submeter_liuwang2023.sh sem jobs liuC na fila. Parar sem perder nada: touch
# $SAIDA/PARAR (as tarefas em curso terminam).
# Saidas por caso em $SAIDA/<id>/: crpc/ (crpc_<a>_<a>_crpc_*.csv e logs) ou mc/ (mc_<a>_<b>.csv, mc_chuva_<a>_<b>.csv);
# eventos em $SAIDA/eventos.log. Resumo: resumo_liuwang2023.py (crpc_liuwang.py).
# Uso: ./submeter_liuwang2023.sh [nos]. Ambiente: SAIDA ($DATA/liuwang2023), M (200), BASE_NYS
# ($DATA/liuwang2023/base_lw_M$M.txt), TRANS (tri3), DT (0.1), PARALELO (24), BLOCO_MC (10), CASOS (todos, ou ids
# separados por virgulas), EXE (bin_liu/SlopeRainMonteCarlo), PYV ($HOME/venv_la/bin/python), INTENS (as 9 do artigo),
# NMC_MAX (sem teto: o N de cada caso mc).
#SBATCH --job-name=liuC
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
M="${M:-200}"
SAIDA="${SAIDA:-$DADOS/liuwang2023}"
BASE_NYS="${BASE_NYS:-$DADOS/liuwang2023/base_lw_M$M.txt}"
TRANS="${TRANS:-tri3}" DT="${DT:-0.1}"
PARALELO="${PARALELO:-${SLURM_CPUS_PER_TASK:-24}}" BLOCO_MC="${BLOCO_MC:-10}"
INTENS="${INTENS:-100,36,15,9,6,4,3,2,1}"
CASOS_TXT="${CASOS_TXT:-$PWD/casos_liuwang2023.txt}" CASOS="${CASOS:-todos}"
EXE="${EXE:-$PWD/bin_liu/SlopeRainMonteCarlo}" PY="$(cd ../../Cho2010 && pwd)" PYV="${PYV:-$HOME/venv_la/bin/python}"
TALUDE=0,0,74,0,74,20,43,20,31,30,0,30
# o caso C: geometria (Fig. 5), nivel inicial reto de (0; 20) a (74; 10) (o chute do permanente) com as cargas de 20 e 10
# m, van Genuchten-Mualem (theta_s = 0,358, theta_r = 0, a = 6,711 kPa: alfa = 0,149 1/kPa, n = 1,289), k_s = 3e-6 m/s,
# q0 = 1200/8760 mm/h, chi por ponto com phi^b = 15 graus, gamma = 20, c' = 5 kPa e phi' = 32 graus (Tabelas 3 e 4),
# chuva por metro inclinado (o fluxo unitario do GeoStudio, setas normais a face na Fig. 5), malhas. Os CoV vao nas opcoes
BASE="geometria=$TALUDE na=14.18918919 na_beta=7.696051722 carga_esq=20 carga_dir=10 vg=0.358,0,0.149011,1.289"
BASE+=" ks=3e-6 q0=0.1369863 chi=ponto phib=15 c=5 phi=32 gamma=20 rcphi=0 transiente=$TRANS dt=$DT fator=tri1 nref=4"
JOB="${SLURM_JOB_ID:-local$$}"
[ -x "$EXE" ] || { echo "falta $EXE: rode ./compilar_liu.sh no no de login" >&2; exit 1; }
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PZ_THREADS=0
mkdir -p "$SAIDA"
LOG="$SAIDA/eventos.log"
diz() { echo "$(date '+%F %T') [$JOB] $*" | tee -a "$LOG"; }

# os casos: id, tipo (crpc|mc), P_F1 do artigo, N (realizacoes; 0: o solo deterministico) e opcoes
IDS=() TIPOS=() NS=() OPCS=()
while read -r id tipo pf n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [ "$CASOS" != todos ] && [[ ",$CASOS," != *",$id,"* ]]; then continue; fi
    [ "$tipo" = mc ] && [ -n "${NMC_MAX:-}" ] && [ "$n" -gt "$NMC_MAX" ] && n="$NMC_MAX"
    IDS+=("$id") TIPOS+=("$tipo") NS+=("$n") OPCS+=("$opc")
done < "$CASOS_TXT"
[ ${#IDS[@]} -gt 0 ] || { echo "nenhum caso (CASOS=$CASOS em $CASOS_TXT)" >&2; exit 1; }

pega_trava() { # como no hu2024_casos.sbatch
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
base_feita() { [ -s "$BASE_NYS" ]; }
NYS_LOCAL=""
base_local() {
    local loc="${TMPDIR:-/tmp}/liuwang2023_$(basename "$BASE_NYS")"
    if [ ! -s "$loc" ] || [ "$(stat -c %s "$loc")" != "$(stat -c %s "$BASE_NYS")" ]; then
        cp "$BASE_NYS" "$loc.$JOB" 2> /dev/null && mv -f "$loc.$JOB" "$loc" || loc="$BASE_NYS"
    fi
    NYS_LOCAL="$loc"
}

# as tarefas: "<i> <rotulo> <opcoes da tarefa>"; feitas quando o log tem 'arquivos:'
tarefas() {
    local i n a b k I
    for ((k = 0; ; k++)); do # CRPC intercaladas: a realizacao k de cada caso
        local alguma=0
        for i in "${!IDS[@]}"; do
            [ "${TIPOS[$i]}" = crpc ] || continue
            n="${NS[$i]}"
            if [ "$n" = 0 ]; then # o solo deterministico: uma tarefa por intensidade, na primeira volta
                if [ "$k" = 0 ]; then
                    for I in ${INTENS//,/ }; do echo "$i det_I$I det=1 n=0 crpc=$I"; done
                    alguma=1
                fi
                continue
            fi
            [ "$k" -lt "$n" ] || continue
            echo "$i r$k n=1 primeira=$k crpc=$INTENS"
            alguma=1
        done
        [ "$alguma" = 1 ] || break
    done
    local maior=0
    for i in "${!IDS[@]}"; do [ "${TIPOS[$i]}" = mc ] && [ "${NS[$i]}" -gt "$maior" ] && maior="${NS[$i]}"; done
    for ((a = 0; a < maior; a += BLOCO_MC)); do # os blocos do MC intercalados entre os casos
        for i in "${!IDS[@]}"; do
            [ "${TIPOS[$i]}" = mc ] && [ "$a" -lt "${NS[$i]}" ] || continue
            b=$((a + BLOCO_MC - 1)); [ "$b" -ge "${NS[$i]}" ] && b=$((NS[i] - 1))
            echo "$i mc_${a}_$b n=$((b - a + 1)) primeira=$a chuva_aleatoria=1"
        done
    done
}
pasta() { [ "${TIPOS[$1]}" = crpc ] && echo "$SAIDA/${IDS[$1]}/crpc" || echo "$SAIDA/${IDS[$1]}/mc"; }
feita() { [ -s "$1/log_$2.txt" ] && grep -q 'arquivos:' "$1/log_$2.txt"; }

# uma onda: ate PARALELO tarefas livres; 1 se nao havia tarefa livre
onda() {
    local n=0 i rot opc d
    base_feita || return 1
    base_local
    while read -r i rot opc; do
        [ "$n" -ge "$PARALELO" ] && break
        [ -e "$SAIDA/PARAR" ] && break
        d="$(pasta "$i")"
        feita "$d" "$rot" && continue
        pega_trava "$d/.tarefas/$rot" || continue
        (
            cd "$d" || exit 1
            # shellcheck disable=SC2086
            if "$EXE" $BASE ${OPCS[$i]} $opc M="$M" nystrom="$NYS_LOCAL" saida="$rot" > "log_$rot.txt" 2>&1; then
                [[ "$rot" == det_* || "$rot" == r* ]] && diz "${IDS[$i]} $rot: $(grep -E 'D_C|FS\(|descartada|parou' "log_$rot.txt" | tr -s ' ' | tr '\n' ';' | cut -c1-300)"
            else
                echo "ERRO $(date '+%F %T')" >> ".tarefas/$rot"
                diz "${IDS[$i]} $rot FALHOU (veja $d/log_$rot.txt)"
            fi
        ) &
        n=$((n + 1))
    done < <(tarefas)
    [ "$n" -gt 0 ] || return 1
    wait
    return 0
}

diz "inicio em $(hostname -s): ${#IDS[@]} casos, M $M, transiente $TRANS (dt $DT h), intensidades $INTENS, PARALELO" \
    "$PARALELO, blocos do MC de $BLOCO_MC, executavel $(cat "$(dirname "$EXE")/versao.txt" 2> /dev/null)"
if ! base_feita && pega_trava "$BASE_NYS.trava"; then
    diz "base de Nystrom $BASE_NYS (M = $M)"
    mkdir -p "$(dirname "$BASE_NYS")"
    if "$PYV" -I "$PY/nystrom_base.py" "$BASE_NYS.tmp" M="$M" talude="$TALUDE" nucleo=exp lh=20 lv=2 dx=0.5 dy=0.1 \
        > "$BASE_NYS.log" 2>&1 && mv -f "$BASE_NYS.tmp" "$BASE_NYS"; then
        diz "base pronta: $(grep -m 1 "nos de" "$BASE_NYS.log" | cut -d, -f1); $(grep -m 1 "M = *$M:" "$BASE_NYS.log" | tr -s ' ')"
    else
        echo "ERRO $(date '+%F %T')" >> "$BASE_NYS.trava"
        diz "base de Nystrom FALHOU (veja $BASE_NYS.log)"
        exit 1
    fi
fi
while :; do
    [ -e "$SAIDA/PARAR" ] && { diz "achei $SAIDA/PARAR: saindo"; break; }
    if ! base_feita; then
        grep -q '^ERRO' "$BASE_NYS.trava" 2> /dev/null && { diz "a base falhou: saindo"; break; }
        sleep 60
        continue
    fi
    onda && continue
    diz "nada mais livre: saindo"
    break
done
diz "fim"
