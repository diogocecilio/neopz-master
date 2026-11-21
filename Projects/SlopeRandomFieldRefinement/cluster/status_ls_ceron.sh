#!/bin/bash
# status_ls_ceron.sh: os jobs lsCeronC na fila, os ultimos eventos e, por caso, o deterministico, o FORM, o line
# sampling (parcial pelo linhas.csv enquanto roda) e o Monte Carlo, contra o artigo (resumo_ls_ceron.py; so python3
# sem numpy: roda no fermi). SAIDA como no submeter_ls_ceron.sh.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/ceron_ls_casos}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n lsCeronC
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_ls_ceron.py "$SAIDA" casos_ls_ceron.txt
