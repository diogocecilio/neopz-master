#!/bin/bash
# baixar_liu2024.sh (no PC): baixa a campanha liu2024_casos.sbatch ($DATA/liu2024 no fermi; pede a senha uma vez) para
# ~/projects/runs/fermi_liu2024 e grava o resumo contra o artigo (resultado.txt, resumo.json; resumo_liu2024.py).
# SEM_BAIXAR=1: so o resumo do que ja foi baixado.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_liu2024}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/liu2024/}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    rsync -avz --exclude '*.xi' --exclude 'l*_*_0_0.csv' --exclude 'it*_mapa_*' --exclude 'mc_mapa_*' "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
python3 -I "$AQUI/resumo_liu2024.py" "$DEST" "$AQUI/casos_liu2024.txt" --json "$DEST/resumo.json" | tee "$DEST/resultado.txt"
