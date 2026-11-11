#!/bin/bash
# status.sh: os jobs slopeMC na fila, a ultima linha do log de cada processo e a tabela do Monte Carlo ate agora (junta os
# mc_resumo.csv e compara cada caso com o artigo). SAIDA como em submeter.sh.
cd "$(dirname "${BASH_SOURCE[0]}")"
SAIDA="${SAIDA:-${DATA:+$DATA/ceron_mc}}"
SAIDA="${SAIDA:-$PWD/resultados}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n slopeMC
if [ -d "$SAIDA/logs" ]; then
    echo "processos: $(ls "$SAIDA"/logs/p*.txt 2>/dev/null | wc -l); tarefas feitas: $(cat "$SAIDA"/logs/p*.txt 2>/dev/null | grep -c 'Gamma =')"
    for f in "$SAIDA"/logs/p*.txt; do [ -e "$f" ] && echo "  $(basename "$f" .txt): $(tail -n 1 "$f")"; done | sort -V | tail -n 5
fi
bin/SlopeMonteCarloCluster resumo "$SAIDA"
