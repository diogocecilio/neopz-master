#!/bin/bash
# submeter_liuwang2023.sh [nos = 4]: submete <nos> jobs liuwang2023_casos.sbatch (um no inteiro cada), que dividem os casos
# de casos_liuwang2023.txt (caso C, Liu e Wang 2023: base de Nystrom, CRPC por realizacao e Monte Carlo direto). Volta em
# segundos: pode desconectar. Pode rodar de novo a qualquer momento: os jobs novos pegam o que estiver livre.
# Sem jobs liuC na fila, libera antes as travas de tarefas que nao terminaram (de jobs que morreram); as com ERRO ficam,
# a menos que RETENTAR=1.
# Ambiente: SAIDA ($DATA/liuwang2023), TEMPO (2-00:00:00), PARTICAO, DEPENDE (id de job: --dependency=afterany: os jobs
# esperam esse terminar), e as do sbatch (M, BASE_NYS, TRANS, DT, PARALELO, BLOCO_MC, INTENS, CASOS), que vao junto
# (--export=ALL). O fermi roda 4 jobs por usuario e aceita 8 na fila: com outras campanhas, menos nos ou DEPENDE.
# FORCAR=1: submete mesmo com o transiente do caso C ainda instavel.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# 22/11: o transiente deste solo (van Genuchten n = 1,289) ainda para em muitas realizacoes (o Newton entra num ciclo
# de 2 iteracoes; no piloto do PC, 1 de 2 realizacoes do MC parou em 2 a 3 h): a campanha so roda com FORCAR=1
if [ -z "${FORCAR:-}" ]; then
    echo "o transiente do caso C ainda nao converge neste solo (veja o README, caso C): FORCAR=1 submete assim mesmo" >&2
    exit 1
fi
NOS="${1:-4}"
case "$NOS" in '' | *[!0-9]* | 0) echo "uso: $0 [numero de nos, padrao 4]" >&2; exit 2 ;; esac
export SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/liuwang2023}"
EXE=bin_liu/SlopeRainMonteCarlo
[ -x "$EXE" ] || { echo "compile antes, no no de login: ./compilar_liu.sh" >&2; exit 1; }
grep -q "crpc_dmax" "$EXE" || { echo "bin_liu/ nao tem as opcoes do caso C: rode ./compilar_liu.sh" >&2; exit 1; }
[ -e "$SAIDA/PARAR" ] && { echo "existe $SAIDA/PARAR: apague antes (rm $SAIDA/PARAR)" >&2; exit 1; }
[ -x "${PYV:-$HOME/venv_la/bin/python}" ] || { echo "falta o python com numpy e scipy (${PYV:-$HOME/venv_la/bin/python}): a base de Nystrom e feita com ele" >&2; exit 1; }
mkdir -p logs "$SAIDA"
fila="$(squeue -h -u "${USER:-$(id -un)}" -n liuC -o %i)" || { echo "o squeue nao respondeu: tente de novo" >&2; exit 1; }
if [ -z "$fila" ]; then
    n=0
    for t in "$SAIDA"/*/*/.tarefas/* "${BASE_NYS:-$SAIDA/base_lw_M${M:-200}.txt}.trava" "$SAIDA"/base_*.trava; do
        [ -e "$t" ] || continue
        if grep -q '^ERRO' "$t" && [ -z "${RETENTAR:-}" ]; then
            echo "  com ERRO (fica; RETENTAR=1 refaz): $t"
            continue
        fi
        case "$t" in
        *.trava) [ -s "${t%.trava}" ] && continue ;;
        *) g="$(dirname "$(dirname "$t")")/log_$(basename "$t").txt"; [ -s "$g" ] && grep -q 'arquivos:' "$g" && continue ;;
        esac
        rm -f "$t"
        n=$((n + 1))
    done
    echo "sem jobs liuC na fila: $n travas de tarefas interrompidas liberadas"
fi
for ((j = 1; j <= NOS; j++)); do
    # shellcheck disable=SC2086
    sbatch --export=ALL ${TEMPO:+--time=$TEMPO} ${PARTICAO:+--partition=$PARTICAO} \
        ${DEPENDE:+--dependency=afterany:$DEPENDE} liuwang2023_casos.sbatch
done
echo "saida: $SAIDA  (acompanhar: ./status_liuwang2023.sh)"
