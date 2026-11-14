#!/bin/bash
# baixar.sh (roda no PC, com a VPN da UFRGS ligada): mostra os jobs chuvaMC no fermi, baixa $DATA/chuva_mc e os logs do
# SLURM, mostra a tabela por caso e gera o relatorio PDF (posproc/). Pede a senha do fermi uma vez so.
#   ./baixar.sh [destino]        destino padrao: ~/projects/runs/fermi_chuva
#   SEM_BAIXAR=1 ./baixar.sh     so o relatorio, com os dados ja baixados (sem VPN)
# O relatorio vai para <destino>/relatorio_parcial enquanto ha jobs na fila e para <destino>/relatorio quando acabaram.
# Ambiente: FERMI (padrao diogo@fermi.cesup.ufrgs.br), REMOTO (padrao /mnt/data/diogo/chuva_mc).
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-$HOME/projects/runs/fermi_chuva}"
HOST="${FERMI:-diogo@fermi.cesup.ufrgs.br}"
REMOTO="${REMOTO:-/mnt/data/diogo/chuva_mc}"
POS="$AQUI/../posproc"
mkdir -p "$DEST"

FILA=""
if [ -z "${SEM_BAIXAR:-}" ]; then
    # uma conexao mestre: a senha e pedida aqui e o ssh e o rsync seguintes a reaproveitam
    TMP="$(mktemp -d)"
    SSH=(ssh -o ControlPath="$TMP/s")
    trap '"${SSH[@]}" -O exit "$HOST" 2> /dev/null; rm -rf "$TMP"' EXIT
    ssh -o ControlMaster=yes -o ControlPath="$TMP/s" -o ControlPersist=yes -fN "$HOST"

    echo "== jobs no fermi =="
    FILA="$("${SSH[@]}" "$HOST" 'squeue -h -u $USER -n chuvaMC -o "%.8i %.2t %.12M %R"')"
    if [ -n "$FILA" ]; then echo "$FILA"; else echo "nenhum job chuvaMC na fila"; fi
    echo "== tabela no fermi =="
    "${SSH[@]}" "$HOST" "SAIDA=$REMOTO neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh" | tail -n +2 || true

    echo "== baixando para $DEST =="
    rsync -az --info=progress2 -e "${SSH[*]}" --exclude .tarefas/ "$HOST:$REMOTO" \
        ":neopz-master/Projects/SlopeRandomFieldRefinement/cluster/logs" "$DEST/" || {
        c=$? # 24: arquivos sumiram durante a copia; 23: alguns arquivos nao foram lidos (ficam as versoes do PC)
        [ $c -eq 23 ] || [ $c -eq 24 ] || { echo "rsync falhou (codigo $c)" >&2; exit $c; }
        echo "AVISO: rsync com codigo $c (veja as mensagens acima); o relatorio usa o que ha no PC" >&2
    }
fi

if [ -n "${SEM_BAIXAR:-}" ] || [ -n "$FILA" ]; then SAIDA="$DEST/relatorio_parcial"; else SAIDA="$DEST/relatorio"; fi
mkdir -p "$SAIDA"
echo "== estatisticas e relatorio em $SAIDA =="
python3 -I -B "$POS/dados_chuva.py" "$DEST/$(basename "$REMOTO")" "$SAIDA/estatisticas.json" > "$SAIDA.log" 2>&1 &&
    python3 -I -B "$POS/relatorio_chuva.py" "$SAIDA/estatisticas.json" "$SAIDA" >> "$SAIDA.log" 2>&1 ||
    { tail -n 20 "$SAIDA.log"; exit 1; }
grep -v "^aviso" "$SAIDA.log" | head -n 40
echo "pronto: $SAIDA/relatorio.pdf  (figuras em $SAIDA/figs, log em $SAIDA.log)"
