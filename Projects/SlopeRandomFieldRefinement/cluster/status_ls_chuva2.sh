#!/bin/bash
# status_ls_chuva2.sh: os jobs lsChuva2 na fila, os ultimos eventos e, por LS, o deterministico, o FORM, o line sampling
# (parcial pelo linhas.csv enquanto roda) e o Monte Carlo direto da chuva2 no mesmo evento (resumo_ls_chuva2.py; so
# python3 sem numpy: roda no fermi). SAIDA como no submeter_ls_chuva2.sh; SAIDA_MC ($DATA/chuva2: o MC, CAMPANHA=chuva2);
# LINHAS (10): os eventos mostrados.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
SAIDA="${SAIDA:-${DATA:-/mnt/data/${USER:-$(id -un)}}/chuva2_ls}"
SAIDA_MC="${SAIDA_MC:-${DATA:-/mnt/data/${USER:-$(id -un)}}/chuva2}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n lsChuva2
[ -f "$SAIDA/PARAR" ] && echo "ATENCAO: existe $SAIDA/PARAR"
if [ -f "$SAIDA/eventos.log" ]; then
    echo "ultimos eventos:"
    tail -n "${LINHAS:-10}" "$SAIDA/eventos.log" | sed 's/^/  /'
fi
python3 -I resumo_ls_chuva2.py "$SAIDA" --mc "$SAIDA_MC"
