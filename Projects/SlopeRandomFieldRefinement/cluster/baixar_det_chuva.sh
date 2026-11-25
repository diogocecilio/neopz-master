#!/bin/bash
# baixar_det_chuva.sh (no PC): baixa a validacao deterministica da chuva (det_chuva.sbatch, $DATA/det_chuva no fermi;
# pede a senha uma vez) para ~/projects/runs/fermi_det_chuva. Pode rodar com os jobs em curso (baixa o que ja existe).
set -euo pipefail
DEST="${DEST:-$HOME/projects/runs/fermi_det_chuva}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/det_chuva/}"
mkdir -p "$DEST"
rsync -avz --exclude '*.tmp' "$ORIGEM" "$DEST/" || { rc=$?; [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"; }
n=$(find "$DEST" -maxdepth 2 -name fim | wc -l)
echo "baixado em $DEST: $n casos terminados"
