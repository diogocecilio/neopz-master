#!/bin/bash
# status.sh: os jobs chuvaMC na fila, as ultimas linhas dos logs dos processos e, por caso, as tarefas feitas, as
# realizacoes ja calculadas e o P_f contado em 24 h ate agora (so com awk: roda no fermi). SAIDA como em submeter.sh.
cd "$(dirname "${BASH_SOURCE[0]}")"
SAIDA="${SAIDA:-${DATA:+$DATA/chuva_mc}}"
SAIDA="${SAIDA:-$PWD/resultados}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n chuvaMC
[ -f "$SAIDA/campanha.txt" ] || { echo "nada em $SAIDA ainda"; exit 0; }
read -r IREF BLOCO < "$SAIDA/campanha.txt"
echo "campanha em $SAIDA: I_ref = $IREF mm/h, blocos de $BLOCO"
if [ -d "$SAIDA/logs" ]; then
    echo "ultimas linhas dos logs:"
    for f in "$SAIDA"/logs/p*.txt; do [ -e "$f" ] && echo "  $(basename "$f" .txt): $(tail -n 1 "$f")"; done | sort -V | tail -n 6
    echo "erros: $(cat "$SAIDA"/logs/p*.txt 2> /dev/null | grep -c ERRO)"
fi
printf "\n%-12s %8s %8s %8s %8s %8s\n" caso tarefas rodando real_24h falhas "P_f(%)"
while read -r id fase n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    d="$SAIDA/$id"
    [ -d "$d" ] || continue
    feitas=$(ls "$d"/mc_mapa_*.csv 2> /dev/null | wc -l)
    travas=$(ls "$d"/.tarefas/ 2> /dev/null | wc -l)
    # realizacoes aleatorias (real >= 0) em 24 h com situacao ok; falha: FS_min < 1
    ls "$d"/mc_[0-9]*_[0-9]*.csv > /dev/null 2>&1 || { printf "%-12s %8d %8d\n" "$id" "$feitas" $((travas - feitas)); continue; }
    awk -F, -v id="$id" -v feitas="$feitas" -v rodando=$((travas - feitas)) '
        FNR == 1 { next }
        $1 >= 0 && $2 == 24 && $4 == "ok" { n++; if ($3 < 1) f++ }
        END { printf "%-12s %8d %8d %8d %8d %8.2f\n", id, feitas, rodando, n, f, n ? 100 * f / n : 0 }' \
        "$d"/mc_[0-9]*_[0-9]*.csv
done < casos.txt
