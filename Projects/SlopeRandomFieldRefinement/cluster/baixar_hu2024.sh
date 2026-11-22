#!/bin/bash
# baixar_hu2024.sh (no PC): baixa a campanha hu2024_casos.sbatch ($DATA/hu2024_ks$KS no fermi; pede a senha uma vez) para
# ~/projects/runs/fermi_hu2024_ks$KS e grava o resumo contra o artigo (resultado.txt, resultado.csv, resumo.json;
# resumo_hu2024.py). SEM_BAIXAR=1: so o resumo do que ja foi baixado. KS (2e-5) escolhe a campanha.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KS="${KS:-2e-5}"
DEST="${DEST:-$HOME/projects/runs/fermi_hu2024_ks$KS}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/hu2024_ks$KS/}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    # as analises em curso (.xi, csv temporarios) nao interessam; o rsync devolve 23/24 quando um arquivo some no meio
    rsync -avz --exclude '*.xi' --exclude 'l*_*_0_0.csv' --exclude 'it*_mapa_*' "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
python3 -I "$AQUI/resumo_hu2024.py" "$DEST" "$AQUI/casos_hu2024.txt" --json "$DEST/resumo.json" \
    --csv "$DEST/resultado.csv" | tee "$DEST/resultado.txt"
