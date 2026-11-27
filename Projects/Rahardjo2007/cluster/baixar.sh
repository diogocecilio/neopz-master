#!/bin/bash
# baixar.sh [modelo] (no PC): baixa as rodadas do Rahardjo2007 do fermi (pede a senha uma vez). modelo iii (o padrao):
# $DATA/rahardjo2007 -> ~/projects/runs/rahardjo2007/fermi; i (ou ii): $DATA/rahardjo2007_i -> .../fermi_i. Pode rodar com
# os jobs em curso (as pastas sem o arquivo fim estao incompletas).
set -euo pipefail
M="${1:-iii}"
case "$M" in iii) SUF="" ;; *[!a-z0-9]*) echo "uso: $0 [iii|ii|i|etapa2]" >&2; exit 2 ;; *) SUF="_$M" ;; esac
DEST="${DEST:-$HOME/projects/runs/rahardjo2007/fermi$SUF}"
ORIGEM="${ORIGEM:-diogo@fermi.cesup.ufrgs.br:/mnt/data/diogo/rahardjo2007$SUF/}"
mkdir -p "$DEST"
rsync -avz "$ORIGEM" "$DEST/" || { rc=$?; [ "$rc" = 23 ] || [ "$rc" = 24 ] || exit "$rc"; }
n=$(ls "$DEST"/*/fim 2> /dev/null | wc -l)
echo "baixado em $DEST: $n casos terminados"
