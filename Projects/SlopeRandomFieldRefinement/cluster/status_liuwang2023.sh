#!/bin/bash
# status_liuwang2023.sh: os jobs liuC na fila, os ultimos eventos e o resumo por caso contra o artigo
# (resumo_liuwang2023.py; so python3 sem numpy: roda no fermi). SAIDA como no submeter_liuwang2023.sh.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/liuwang2023}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n liuC
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_liuwang2023.py "$SAIDA" casos_liuwang2023.txt
