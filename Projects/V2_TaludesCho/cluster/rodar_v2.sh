#!/bin/bash
# rodar_v2.sh <primeira> <ultima>: roda as realizacoes primeira..ultima do Monte Carlo do V2_TaludesCho, uma execucao do
# programa por realizacao (o contrato: REAL, CAMPOS, CASO, NUCLEO, ...), cada uma em <SAIDA>/<COMB>/r<nnnnn>/, com
# PARALELO execucoes ao mesmo tempo (padrao: nproc; cada uma usa 1 CPU). Retomavel: as que ja tem resultado.txt sao
# puladas (resultado.txt e gravado por ultimo). Uma realizacao que sai com codigo != 0 vai para <COMB>/FALHAS e e
# tentada de novo na proxima chamada. DET=1 roda antes o caso deterministico (CAMPOS=0) em <COMB>/det/.
# Ambiente: EXE (padrao ./bin/V2_TaludesCho), SAIDA (padrao ./resultados), CASO (cho1|cho2), NUCLEO (gauss|exp), COMB
# (nome da pasta; padrao <CASO>_<NUCLEO>) e qualquer troca do executavel (COV_C, COV_PHI, LX, LY, RHO, KL_FRAC,
# H_MALHA, SEMENTE, FS_MAX, ...): use um COMB proprio para cada conjunto de trocas. No fim: ./pf_v2.py <SAIDA>/<COMB>.
set -u
export LC_ALL=C
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -eq 2 ] || { sed -n '2,/^set -u/{/^#/p}' "$0"; exit 2; }
export EXE="${EXE:-$PKG/bin/V2_TaludesCho}"
[ -x "$EXE" ] || { echo "executavel $EXE nao encontrado (rode ./compilar_v2.sh)" >&2; exit 2; }
export CASO="${CASO:-cho1}" NUCLEO="${NUCLEO:-gauss}" MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 VERBOSO="${VERBOSO:-0}"
export DIR_MC="${SAIDA:-$PKG/resultados}/${COMB:-${CASO}_${NUCLEO}}"
export V2_VERSAO="${V2_VERSAO:-$(cat "$(dirname "$EXE")/versao.txt" 2>/dev/null || echo -)}"
mkdir -p "$DIR_MC" || exit 1
env | grep -E '^(CASO|NUCLEO|SEMENTE|COV_C|COV_PHI|LX|LY|RHO|KL_FRAC|KL_M|KL_MARGEM|H_MALHA|ORDEM|FS_TOL|FS_F0|FS_MAX|MAXIT_FROZEN|COES|PHI_DEG|E|NU|GAMMA|V2_VERSAO)=' \
    | sort > "$DIR_MC/parametros.txt"
uma() { # uma realizacao (r = numero) ou o deterministico (r = det)
    local r="$1" d e
    if [ "$r" = det ]; then d="$DIR_MC/det"; else d="$DIR_MC/$(printf 'r%05d' "$r")"; fi
    [ -f "$d/resultado.txt" ] && return 0
    mkdir -p "$d" && cd "$d" || return 0
    if [ "$r" = det ]; then CAMPOS=0 REAL=0 "$EXE" > log.txt 2>&1; else CAMPOS=1 REAL="$r" "$EXE" > log.txt 2>&1; fi
    e=$?
    [ $e -eq 0 ] || echo "$(date '+%F %T') $(hostname -s) realizacao $r: codigo de saida $e" >> "$DIR_MC/FALHAS"
    return 0
}
export -f uma
[ "${DET:-0}" = 1 ] && uma det
seq "$1" "$2" | xargs -P "${PARALELO:-$(nproc)}" -I{} bash -c 'uma {}'
echo "$(date '+%F %T') $DIR_MC: $(ls "$DIR_MC"/r*/resultado.txt 2>/dev/null | wc -l) realizacoes concluidas"
