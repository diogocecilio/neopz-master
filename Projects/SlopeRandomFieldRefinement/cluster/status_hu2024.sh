#!/bin/bash
# status_hu2024.sh: os jobs huA na fila, os ultimos eventos e, por caso, o deterministico, o FORM, o line sampling (parcial
# pelo linhas.csv enquanto roda) e o Monte Carlo, contra o artigo (resumo_hu2024.py; so python3 sem numpy: roda no
# fermi). KS e SAIDA como no submeter_hu2024.sh.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
KS="${KS:-2e-5}"
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/hu2024_ks$KS}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n huA
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_hu2024.py "$SAIDA" casos_hu2024.txt
