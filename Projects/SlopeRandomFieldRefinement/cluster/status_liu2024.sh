#!/bin/bash
# status_liu2024.sh: os jobs liuB na fila, os ultimos eventos e, por cenario, o deterministico, o FORM, o line sampling
# (parcial enquanto roda) e o Monte Carlo, contra o artigo (resumo_liu2024.py; so python3 sem numpy: roda no fermi).
# SAIDA como no submeter_liu2024.sh.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/liu2024}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n liuB
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_liu2024.py "$SAIDA" casos_liu2024.txt
