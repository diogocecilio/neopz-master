#!/bin/bash
# submeter_mc.sh <projeto> [nos]: o Monte Carlo de um projeto (V2_TaludesCho ou V4_TaludeLiu; V2 e V4 bastam; ou
# PROJETO=...) no SLURM (fermi). Compila, se preciso, no no de login em segundo plano (../compilar.sh), e submete <nos>
# jobs job_mc.sbatch (padrao 4), um no inteiro cada, que dividem a mesma fila: o caso deterministico de cada entrada
# (DET=1, o padrao) e depois as realizacoes PRIMEIRA..ULTIMA (padrao: as do artigo, em projetos.sh) de todas as entradas
# de CASOS (padrao: a campanha do artigo), intercaladas. Volta em segundos: pode desconectar.
# Rodar de novo submete mais nos para a mesma fila (o que esta concluido ou rodando nao se repete), e e tambem o jeito
# de retomar uma campanha parada (apaga $SAIDA/<projeto>/PARAR). Um job que chega perto do limite de tempo (TEMPO,
# padrao 4-00:00:00) submete a sua continuacao. Se o SLURM recusar um job (cota de jobs por usuario), os ja submetidos
# seguem. Os avisos (cota, parametros, bibliotecas, falhas) ficam em $SAIDA/<projeto>/AVISOS, que status_mc.sh mostra.
# Ambiente: CASOS, PRIMEIRA, ULTIMA, DET, SAIDA (padrao $DATA/monte_carlo), TEMPO, PARTICAO, RECOMPILAR=1, FORCAR=1
# (recompila mesmo com jobs na fila), FORCAR_SEM_COMPILAR=1 (submete com o executavel atual, sem recompilar) e as
# variaveis de rodar_mc.sh e do programa (COV_C, H_MALHA, ...: valem para todas as entradas), que vao para os jobs.
set -euo pipefail
MC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$MC/.." && pwd)"
cd "$PKG"
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
source "$MC/projetos.sh"
if [ $# -gt 0 ] && ! [[ "$1" =~ ^[0-9]+$ ]]; then PROJETO="$1"; shift; fi
[ -n "${PROJETO:-}" ] || { echo "uso: $0 <projeto> [numero de nos, padrao 4] (projetos: $(mc_projetos))" >&2; exit 2; }
mc_projeto "$PROJETO" || exit 2
export PROJETO="$MC_PROJETO"
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 <projeto> [numero de nos, padrao 4]" >&2; exit 2 ;; esac
recusa_se_compilando
# resultados: $DATA/monte_carlo/<projeto>, na particao de dados de onde o manual do fermi manda rodar os jobs
# (temporaria, sem backup: copie para o PC ou para o $HOME ao terminar); sem $DATA, resultados_mc nesta pasta
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
export SAIDA="${SAIDA:-$PKG/resultados_mc}" VERIF_CLUSTER="$PKG"
CASOS="$(echo ${CASOS:-$MC_CASOS_PADRAO})"
export CASOS PRIMEIRA="${PRIMEIRA:-$MC_PRIMEIRA}" ULTIMA="${ULTIMA:-$MC_ULTIMA}" DET="${DET:-1}"
[ -n "${TEMPO:-}" ] && export TEMPO
[[ "$PRIMEIRA" =~ ^[0-9]+$ && "$ULTIMA" =~ ^[0-9]+$ && "$PRIMEIRA" -le "$ULTIMA" ]] ||
    { echo "PRIMEIRA=$PRIMEIRA e ULTIMA=$ULTIMA: inteiros com PRIMEIRA <= ULTIMA" >&2; exit 2; }
# a mesma validacao de rodar_mc.sh (o modo --pastas so le CASOS), antes de compilar e submeter
"$MC/rodar_mc.sh" --pastas > /dev/null || exit 2
D="$SAIDA/$MC_PROJETO"
mkdir -p "$PKG/logs" "$D"
JOB="verif-mc-$MC_ID"
aviso() { echo "$(date '+%F %T') $*" | tee -a "$D/AVISOS" >&2; }
# uma variavel que o programa le e projetos.sh nao conhece ficaria fora de parametros.txt
novas="$(mc_confere_vars "$VERIF_NEOPZ/Projects/$MC_PROJETO" | tr '\n' ' ')"
[ -z "$novas" ] || echo "aviso: $MC_PROJETO le variaveis que mc/projetos.sh nao classifica: $novas(nao entram em" \
    "parametros.txt; acrescente-as la em MC_VARS ou MC_VARS_CONTROLE)" >&2

# a campanha (para status_mc.sh), ja antes da compilacao
{
    echo "PROJETO=$MC_PROJETO"
    echo "CASOS=\"$CASOS\""
    echo "PRIMEIRA=$PRIMEIRA"
    echo "ULTIMA=$ULTIMA"
    echo "DET=$DET"
} > "$D/campanha.txt"

# PARAR: submeter e um pedido de recomeco. O de todos os projetos ($SAIDA/PARAR) so sai se nenhum outro projeto tem jobs
# na fila (eles voltariam a rodar); a cadeia da compilacao respeita um PARAR criado durante ela (o usuario parou depois)
tira_parar() {
    if [ -e "$SAIDA/PARAR" ]; then
        local outros por=""
        if ! outros="$(jobs_verif 2> /dev/null | awk -v n="$JOB" '$2 ~ /^verif-mc-/ && $2 != n' | wc -l)"; then
            outros="?" por=" e o squeue nao respondeu (pode haver jobs de outros projetos na fila)"
        elif [ "$outros" -gt 0 ]; then
            por=" e ha $outros job(s) de outros projetos na fila (voltariam a rodar)"
        fi
        if [ -n "${CADEIA:-}" ] || [ "$outros" != 0 ]; then
            aviso "$SAIDA/PARAR (a parada de todos os projetos) existe$por: nenhum job submetido. Apague-o" \
                "(rm $SAIDA/PARAR) e rode mc/submeter_mc.sh $MC_ID $NOS"
            exit 1
        fi
        rm -f "$SAIDA/PARAR"
        echo "apagado $SAIDA/PARAR (a parada de todos os projetos)"
    fi
    if [ -e "$D/PARAR" ]; then
        if [ -n "${CADEIA:-}" ]; then
            aviso "compilacao terminada, mas $D/PARAR existe: nenhum job submetido (apague-o e rode mc/submeter_mc.sh $MC_ID $NOS)"
            exit 1
        fi
        rm -f "$D/PARAR"
        echo "apagado $D/PARAR: a campanha volta a rodar"
    fi
}
# este pedido e um pedido de recomeco: um PARAR de antes sai ja, tambem quando e preciso compilar antes
[ -z "${CADEIA:-}" ] && tira_parar
# compilacao no no de login, em segundo plano, se o executavel nao existe ou e mais velho que os fontes; ao fim ela
# repete este pedido (CADEIA=1), que entao submete
export AVISOS_CADEIA="$D/AVISOS"
compilar_antes "$MC_PROJETO" "$NOS no(s), $MC_PROJETO, CASOS=\"$CASOS\", realizacoes $PRIMEIRA..$ULTIMA, DET=$DET, SAIDA=$SAIDA" \
    "$MC/submeter_mc.sh" "$NOS"
[ -n "${CADEIA:-}" ] && tira_parar

# os parametros de cada pasta conferem com os de parametros.txt (senao os jobs sairiam logo no inicio)
if ! verif="$("$MC/rodar_mc.sh" --verifica 0 0 2>&1 | tr '\n' ' '; exit "${PIPESTATUS[0]}")"; then
    aviso "nao submetido: $verif"
    exit 2
fi
# travas de parada de jobs que ja acabaram
for f in "$D"/.parar.*; do
    [ -e "$f" ] || continue
    j="${f##*.}"
    if out="$(squeue -h -j "$j" -o %i 2>&1)"; then [ -n "$out" ] || rm -f "$f"
    else case "$out" in *[Ii]"nvalid job id"*) rm -f "$f" ;; esac
    fi
done
ja="$(jobs_verif 2> /dev/null | awk -v n="$JOB" '$2 == n' | wc -l)" || ja=0
[ "$ja" -gt 0 ] && echo "ja ha $ja job(s) $JOB na fila"
P=()
[ -n "${PARTICAO:-}" ] && P=(--partition="$PARTICAO")
ok=0
err="$(mktemp)"
for ((i = 1; i <= NOS; i++)); do
    if out="$(sbatch --parsable --job-name="$JOB" --output="$PKG/logs/%x_%j.out" ${P[@]+"${P[@]}"} \
        ${TEMPO:+--time="$TEMPO"} "$MC/job_mc.sbatch" 2> "$err")"; then
        j="${out##*$'\n'}"
        j="${j%%;*}"
        ok=$((ok + 1))
        echo "job $j submetido"
        echo "$(date '+%F %T') $j $MC_PROJETO CASOS=\"$CASOS\" realizacoes $PRIMEIRA..$ULTIMA DET=$DET SAIDA=$SAIDA" \
            >> "$PKG/logs/jobs_mc.txt"
    else
        aviso "o SLURM recusou o job $i de $NOS: $(cat "$err") $out"
        break
    fi
done
rm -f "$err"
[ "$ok" -gt 0 ] && rm -f "$D/CONTINUACAO_PENDENTE"
echo "$ok de $NOS jobs submetidos: $MC_PROJETO, $(echo $CASOS | wc -w) entrada(s), realizacoes $PRIMEIRA..$ULTIMA$(
    [ "$DET" = 1 ] && echo ' e os deterministicos'), resultados em $D"
echo "acompanhe com: mc/status_mc.sh $MC_ID"
[ "$ok" -gt 0 ]
