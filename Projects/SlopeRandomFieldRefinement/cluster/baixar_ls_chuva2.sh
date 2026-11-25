#!/bin/bash
# baixar_ls_chuva2.sh (no PC): baixa o FORM + line sampling da chuva2 (ls_chuva2.sbatch, $DATA/chuva2_ls no fermi; pede a
# senha uma vez) para ~/projects/runs/fermi_chuva2_ls e grava o resumo contra o Monte Carlo direto da chuva2 no mesmo
# evento (resultado.txt, resumo.json; resumo_ls_chuva2.py). O MC vem do que o baixar.sh (CAMPANHA=chuva2) ja baixou
# (MC, padrao ~/projects/runs/fermi_chuva2/chuva2; sem ele, o resumo sai sem o MC). SEM_BAIXAR=1: so o resumo do que ja
# foi baixado. Sem os .xi, os csv das analises em curso e os arquivos do transiente das linhas em voo.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_chuva2_ls}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/chuva2_ls/}"
MC="${MC:-$HOME/projects/runs/fermi_chuva2/chuva2}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    # o rsync devolve 23/24 quando um arquivo some no meio (as analises em curso)
    rsync -avz --exclude '*.xi' --exclude 'l*_*_0_0.csv' --exclude 'l*_*_parametros.txt' --exclude 'l*_pw.txt*' \
        --exclude 'sonda_reuso*' --exclude 'it*_mapa_*' "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
args=("$DEST" --json "$DEST/resumo.json")
[ -d "$MC" ] && args+=(--mc "$MC")
python3 -I "$AQUI/resumo_ls_chuva2.py" "${args[@]}" | tee "$DEST/resultado.txt"
