#!/bin/bash
# baixar_ls_ceron.sh (no PC): baixa a campanha ls_ceron_casos.sbatch ($DATA/ceron_ls_casos no fermi; pede a senha uma
# vez) para ~/projects/runs/fermi_ceron_ls_casos, grava o resumo contra o artigo (resultado.txt, resultado.csv,
# resumo.json; resumo_ls_ceron.py) e as figuras da Secao 6 do artigo (figs/; figuras_ls_ceron.py, com numpy e
# matplotlib). SEM_BAIXAR=1: so o resumo e as figuras do que ja foi baixado.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_ceron_ls_casos}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/ceron_ls_casos/}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    # as analises em curso (.xi, csv temporarios) nao interessam; o rsync devolve 23/24 quando um arquivo some no meio
    rsync -avz --exclude '*.xi' --exclude 'l*_*_0_0.csv' --exclude 'it*_mapa_*' "$ORIGEM" "$DEST/" || {
        rc=$?
        [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"
    }
fi
python3 -I "$AQUI/resumo_ls_ceron.py" "$DEST" "$AQUI/casos_ls_ceron.txt" --json "$DEST/resumo.json" \
    --csv "$DEST/resultado.csv" | tee "$DEST/resultado.txt"
python3 -I "$(cd "$AQUI/../../SlopeRandomField/comparacao_ceron" && pwd)/figuras_ls_ceron.py" "$DEST/resumo.json" \
    "$DEST/figs" > "$DEST/figs.log" 2>&1 && echo "figuras em $DEST/figs" || { echo "as figuras falharam:"; tail -n 20 "$DEST/figs.log"; }
