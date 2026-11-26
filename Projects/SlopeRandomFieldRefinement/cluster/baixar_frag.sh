#!/bin/bash
# baixar_frag.sh (no PC): baixa a campanha frag ($DATA/frag no fermi; pede a senha uma vez) para
# ~/projects/runs/fermi_frag e roda a verificacao da curva de fragilidade (resumo_frag.py): resultado.txt, resumo.json e
# as figuras em figs/. Pode rodar com os jobs em curso (usa o que ja terminou). SEM_BAIXAR=1: so o resumo.
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${DEST:-$HOME/projects/runs/fermi_frag}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/frag/}"
mkdir -p "$DEST"
if [ -z "${SEM_BAIXAR:-}" ]; then
    rsync -avz --exclude '*_real-1/' --exclude .tarefas/ "$ORIGEM" "$DEST/" || { rc=$?; [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"; }
fi
python3 -I "$AQUI/resumo_frag.py" "$DEST" --json "$DEST/resumo.json" --figs "$DEST/figs" | tee "$DEST/resultado.txt"
