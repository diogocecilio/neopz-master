#!/bin/bash
# baixar_liuwang2023.sh (no PC): baixa a campanha liuwang2023_casos.sbatch ($DATA/liuwang2023 no fermi; pede a senha uma
# vez) para ~/projects/runs/fermi_liuwang2023 e grava o resumo contra o artigo (resultado.txt, resumo.json;
# resumo_liuwang2023.py). SEM_BAIXAR=1: so o resumo do que ja foi baixado.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_liuwang2023}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/liuwang2023/}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    rsync -avz --exclude '*_mapa_*' --exclude 'base_*.txt' "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
python3 -I "$AQUI/resumo_liuwang2023.py" "$DEST" "$AQUI/casos_liuwang2023.txt" --json "$DEST/resumo.json" | tee "$DEST/resultado.txt"
