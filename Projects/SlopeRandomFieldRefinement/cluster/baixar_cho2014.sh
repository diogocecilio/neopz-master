#!/bin/bash
# baixar_cho2014.sh (no PC): baixa a campanha cho2014_casos.sbatch ($DATA/cho2014 no fermi; pede a senha uma vez) para
# ~/projects/runs/fermi_cho2014 e grava o resumo contra o artigo e contra a referencia em Python (resultado.txt,
# resumo.json; resumo_cho2014.py; a referencia em ~/projects/runs/cho2014_python, se existir). SEM_BAIXAR=1: so o
# resumo do que ja foi baixado. Sem as bases de Nystrom (refeitas em segundos) e sem os .xi do FORM.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_cho2014}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/cho2014/}"
PYREF="${PYREF:-$HOME/projects/runs/cho2014_python}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    rsync -avz --exclude '*.xi' --exclude 'bases/*.txt' --exclude 'it*_0_0.csv' --exclude 'it*_parametros.txt' \
        "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
args=("$DEST")
[ -d "$PYREF" ] && args+=(--python "$PYREF")
python3 -I "$AQUI/resumo_cho2014.py" "${args[@]}" --json "$DEST/resumo.json" | tee "$DEST/resultado.txt"
