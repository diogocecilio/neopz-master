#!/bin/bash
# submeter_caso.sh <projeto> <caso>[:VAR=valor,...] [...]: submete um job caso.sbatch (um no inteiro) com os casos
# completos de um projeto (V1..V6, ou o nome completo), que rodam ao mesmo tempo, cada um na sua pasta
# $SAIDA_CASOS/<projeto>/<caso>/ (padrao $DATA/verificacao). O nome do caso e o da pasta e, sem CASO nas trocas, o CASO
# do programa. Compila antes, se preciso, no no de login em segundo plano (como mc/submeter_mc.sh), e submete ao fim.
# Exemplos (os casos de cada projeto estao no README.md dele):
#   ./submeter_caso.sh V1 todos
#   ./submeter_caso.sh V3 retangular inclinada dois_blocos retangular_a2:CASO=retangular,AG=2
#   ./submeter_caso.sh V5 tabela1:VTK=1 tabela4:NU_LISTA=0.2,0.25,0.3,0.35,0.4 tabela4_ks001:CASO=tabela4,KS_INF=0.01
#   ./submeter_caso.sh V6 tabelas b40_h10:N=100000,SEMENTE=101 fig15:N=100000,SEMENTE=110
# Ambiente: SAIDA_CASOS, TEMPO (limite do job, padrao 1-00:00:00), PARTICAO, PARALELO, REFAZER=1 (roda de novo os
# concluidos), RECOMPILAR=1, FORCAR=1, FORCAR_SEM_COMPILAR=1 (ver mc/submeter_mc.sh) e variaveis do programa para todos
# os casos (ficam fora de variaveis.txt: prefira as trocas).
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PKG"
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
[ $# -ge 2 ] || { sed -n '2,/^set -e/{/^#/p}' "$0"; exit 2; }
PJ="$(projeto_completo "$1")" || exit 2
shift
vistos=" "
for e in "$@"; do
    m="$(valida_caso "$e")"
    [ -z "$m" ] || { echo "$m" >&2; exit 2; }
    [[ "$vistos" == *" ${e%%:*} "* ]] && { echo "o caso ${e%%:*} aparece duas vezes (cada caso e uma pasta)" >&2; exit 2; }
    vistos="$vistos${e%%:*} "
done
recusa_se_compilando
SAIDA_CASOS="${SAIDA_CASOS:-${DATA:+$DATA/verificacao}}"
export SAIDA_CASOS="${SAIDA_CASOS:-$PKG/resultados_casos}" VERIF_CLUSTER="$PKG" PROJETO="$PJ" CASOS="$*"
mkdir -p "$PKG/logs"
export AVISOS_CADEIA="$PKG/logs/avisos.txt"
compilar_antes "$PJ" "casos de $PJ: $*" "$PKG/submeter_caso.sh" "$PJ" "$@"
ID="${PJ%%_*}"
A=()
[ -n "${PARTICAO:-}" ] && A+=(--partition="$PARTICAO")
[ -n "${TEMPO:-}" ] && A+=(--time="$TEMPO")
err="$(mktemp)"
if out="$(sbatch --parsable --job-name="verif-caso-$ID" --output="$PKG/logs/%x_%j.out" ${A[@]+"${A[@]}"} \
    "$PKG/caso.sbatch" 2> "$err")"; then
    j="${out##*$'\n'}"
    j="${j%%;*}"
    rm -f "$err"
    echo "$(date '+%F %T') $j $PJ CASOS=\"$*\" SAIDA_CASOS=$SAIDA_CASOS" >> "$PKG/logs/jobs_caso.txt"
    echo "job $j submetido: $PJ, $# caso(s), resultados em $SAIDA_CASOS/$PJ/<caso>/; log em logs/verif-caso-${ID}_$j.out"
    echo "acompanhe com: ./status.sh"
else
    echo "$(date '+%F %T') o SLURM recusou o job de $PJ ($*): $(cat "$err") $out" | tee -a "$PKG/logs/avisos.txt" >&2
    rm -f "$err"
    exit 1
fi
