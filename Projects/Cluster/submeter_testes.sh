#!/bin/bash
# submeter_testes.sh [projeto ...]: submete um job testes.sbatch (um no inteiro): os testes padrao (ctest) e os
# [.lento] dos projetos (padrao: todos), com os logs em logs/testes_<job>/. Compila antes, se preciso, no no de login em
# segundo plano (como mc/submeter_mc.sh), e submete ao fim (se um projeto nao compilou, o job o testa assim mesmo e avisa
# que o teste e mais velho que os fontes).
# Ambiente: LENTOS=0 (so os padrao), TEMPO (limite do job, padrao 12:00:00), PARTICAO, TIMEOUT_TESTE, as variaveis dos
# testes (V2_NMC, ...), RECOMPILAR=1, FORCAR=1, FORCAR_SEM_COMPILAR=1.
set -euo pipefail
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PKG"
source "$PKG/modulos.sh"
source "$PKG/comum.sh"
PS=()
if [ $# -gt 0 ]; then
    for a in "$@"; do p="$(projeto_completo "$a")" || exit 2; PS+=("$p"); done
else
    mapfile -t PS < <(projetos_verif)
fi
recusa_se_compilando
export VERIF_CLUSTER="$PKG"
[ $# -gt 0 ] && export PROJETOS="${PS[*]}"
mkdir -p "$PKG/logs"
export AVISOS_CADEIA="$PKG/logs/avisos.txt" CADEIA_SEGUE=1
compilar_antes "${PS[*]}" "testes de ${PS[*]}" "$PKG/submeter_testes.sh" "$@"
A=()
[ -n "${PARTICAO:-}" ] && A+=(--partition="$PARTICAO")
[ -n "${TEMPO:-}" ] && A+=(--time="$TEMPO")
err="$(mktemp)"
if out="$(sbatch --parsable --job-name=verif-testes --output="$PKG/logs/%x_%j.out" ${A[@]+"${A[@]}"} \
    "$PKG/testes.sbatch" 2> "$err")"; then
    j="${out##*$'\n'}"
    j="${j%%;*}"
    rm -f "$err"
    echo "$(date '+%F %T') $j testes de ${PS[*]} LENTOS=${LENTOS:-1}" >> "$PKG/logs/jobs_testes.txt"
    echo "job $j submetido: testes de ${PS[*]}$([ "${LENTOS:-1}" = 1 ] && echo ' (padrao e [.lento])');" \
        "resumo em logs/testes_$j/resumo.txt"
else
    echo "$(date '+%F %T') o SLURM recusou o job de testes: $(cat "$err") $out" | tee -a "$PKG/logs/avisos.txt" >&2
    rm -f "$err"
    exit 1
fi
