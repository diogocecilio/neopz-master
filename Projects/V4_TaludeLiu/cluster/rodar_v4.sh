#!/bin/bash
# rodar_v4.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo do V4_TaludeLiu, uma execucao do
# programa por realizacao (o contrato: CASO, ANGULO, REAL, CAMPOS, SEMENTE, ...), cada uma em <SAIDA>/<COMB>/r<nnnnn>/,
# com PARALELO execucoes ao mesmo tempo (padrao: nproc; cada uma usa 1 CPU). Retomavel: as que ja tem resultado.txt sao
# puladas (resultado.txt e gravado por ultimo). Uma realizacao que sai com codigo != 0 vai para <COMB>/FALHAS e e
# tentada de novo na proxima chamada. DET=1 roda antes o caso deterministico (CASO=det) em <COMB>/det/. No cenario s2
# (k_s constante) o regime permanente e o do deterministico: o det roda primeiro e as realizacoes leem o seu pw.csv
# (PW_PERMANENTE), o que poupa o transiente.
# Ambiente: EXE (padrao ./bin/V4_TaludeLiu), SAIDA (padrao ./resultados), CASO (s1|s2|s3, padrao s3), ANGULO (0|90),
# COMB (nome da pasta; padrao <CASO>_a<ANGULO>, com _rho<RHO> se RHO for dado) e qualquer troca do executavel (RHO,
# COV_C, COV_PHI, COV_K, LX, LY, KL_M, MALHA, P_KR, Q_KS, SEMENTE, MASSA_LIMIAR, ...): use um COMB proprio para cada
# conjunto de trocas. No fim: ./estatisticas_v4.py <SAIDA>/<COMB>.
set -u
export LC_ALL=C
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -eq 2 ] || { sed -n '2,/^set -u/{/^#/p}' "$0"; exit 2; }
export EXE="${EXE:-$PKG/bin/V4_TaludeLiu}"
[ -x "$EXE" ] || { echo "executavel $EXE nao encontrado (rode ./compilar_v4.sh)" >&2; exit 2; }
export CASO="${CASO:-s3}" ANGULO="${ANGULO:-0}" MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 VERBOSO="${VERBOSO:-0}"
case "$CASO" in s1|s2|s3) ;; *) echo "CASO deve ser s1, s2 ou s3 (o deterministico vai com DET=1)" >&2; exit 2 ;; esac
case "$ANGULO" in 0|90) ;; *) echo "ANGULO deve ser 0 ou 90" >&2; exit 2 ;; esac
COMB_PADRAO="${CASO}_a${ANGULO}"
[ -n "${RHO:-}" ] && COMB_PADRAO="${COMB_PADRAO}_rho${RHO}"
export DIR_MC="${SAIDA:-$PKG/resultados}/${COMB:-$COMB_PADRAO}"
export V4_VERSAO="${V4_VERSAO:-$(cat "$(dirname "$EXE")/versao.txt" 2>/dev/null || echo -)}"
mkdir -p "$DIR_MC" || exit 1
DIR_MC="$(cd "$DIR_MC" && pwd)" # absoluto: cada realizacao roda na sua pasta
env | grep -E '^(CASO|ANGULO|SEMENTE|COV_C|COV_PHI|COV_K|RHO|LX|LY|THETA1|THETA2|KL_M|KL_FRAC|MALHA|P_KR|Q_KS|KS|COES|PHI_DEG|PHIB_DEG|E|NU|GAMMA|FS_TOL|FS_F0|MAXIT_FROZEN|MASSA_LIMIAR|TOL_PERM|V4_VERSAO)=' \
    | sort > "$DIR_MC/parametros.txt"
uma() { # uma realizacao (r = numero) ou o deterministico (r = det)
    local r="$1" d e
    if [ "$r" = det ]; then d="$DIR_MC/det"; else d="$DIR_MC/$(printf 'r%05d' "$r")"; fi
    [ -f "$d/resultado.txt" ] && return 0
    mkdir -p "$d" && cd "$d" || return 0
    if [ "$r" = det ]; then CASO=det CAMPOS=0 REAL=0 PW_PERMANENTE= "$EXE" > log.txt 2>&1
    else CAMPOS=1 REAL="$r" "$EXE" > log.txt 2>&1; fi
    e=$?
    [ $e -eq 0 ] || echo "$(date '+%F %T') $(hostname -s) realizacao $r: codigo de saida $e" >> "$DIR_MC/FALHAS"
    return 0
}
export -f uma
if [ "$CASO" = s2 ] && [ -z "${PW_PERMANENTE:-}" ]; then
    (uma det)
    [ -f "$DIR_MC/det/pw.csv" ] && export PW_PERMANENTE="$DIR_MC/det/pw.csv"
elif [ "${DET:-0}" = 1 ]; then
    (uma det)
fi
seq "$1" "$2" | xargs -P "${PARALELO:-$(nproc)}" -I{} bash -c 'uma {}'
echo "$(date '+%F %T') $DIR_MC: $(ls "$DIR_MC"/r*/resultado.txt 2>/dev/null | wc -l) realizacoes concluidas"
