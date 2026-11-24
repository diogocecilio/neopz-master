#!/bin/bash
# status_cho2014.sh: os jobs cho14 na fila, os ultimos eventos e o resumo da campanha cho2014_casos.sbatch por caso
# (deterministico, Monte Carlo com as P_f de 18 e 24 h contra a Fig. 10, FORM e line sampling; resumo_cho2014.py, so
# python3 sem numpy: roda no fermi). SAIDA como no submeter_cho2014.sh; LINHAS (10): os eventos mostrados.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/cho2014}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n cho14
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_cho2014.py "$SAIDA" --curto
