#!/bin/bash
# baixar.sh (roda no PC, com a VPN da UFRGS ligada): mostra os jobs slopeMC no fermi, baixa a campanha e os logs do
# SLURM, mostra a tabela do Monte Carlo e gera o relatorio PDF. Pede a senha do fermi uma vez so.
#   ./baixar.sh [destino]        a campanha de 4 niveis ($DATA/ceron_mc_h4); destino padrao: ~/projects/runs/fermi_ceron
#   NREF=3 ./baixar.sh           a de 3 niveis ($DATA/ceron_mc)
#   SEM_BAIXAR=1 ./baixar.sh     so a tabela e o relatorio, com os dados ja baixados (sem VPN)
# O relatorio vai para <destino>/relatorio_parcial_h<NREF> enquanto ha jobs na fila e para <destino>/relatorio_h<NREF>
# quando acabaram (na de 3 niveis, relatorio_parcial e relatorio).
set -euo pipefail
NREF="${NREF:-4}"
if [ "$NREF" = 3 ]; then PASTA=ceron_mc SUF=""; else PASTA="ceron_mc_h$NREF" SUF="_h$NREF"; fi
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-$HOME/projects/runs/fermi_ceron}"
HOST="${FERMI:-diogo@fermi.cesup.ufrgs.br}"
BIN="${BIN:-$HOME/projects/neopz-master-build-release/Projects/SlopeMonteCarlo/SlopeMonteCarloCluster}"
SCRIPTS="$AQUI/../../SlopeRandomField/comparacao_ceron"
mkdir -p "$DEST"

FILA=""
if [ -z "${SEM_BAIXAR:-}" ]; then
    # uma conexao mestre: a senha e pedida aqui e o ssh e o rsync seguintes a reaproveitam
    TMP="$(mktemp -d)"
    SSH=(ssh -o ControlPath="$TMP/s")
    trap '"${SSH[@]}" -O exit "$HOST" 2> /dev/null; rm -rf "$TMP"' EXIT
    ssh -o ControlMaster=yes -o ControlPath="$TMP/s" -o ControlPersist=yes -fN "$HOST"

    echo "== jobs no fermi =="
    FILA="$("${SSH[@]}" "$HOST" 'squeue -h -u $USER -n slopeMC -o "%.8i %.2t %.12M %R"')"
    if [ -n "$FILA" ]; then echo "$FILA"; else echo "nenhum job slopeMC na fila: a simulacao terminou"; fi

    echo "== baixando para $DEST =="
    rsync -az --info=progress2 -e "${SSH[*]}" --exclude danificados/ "$HOST:/mnt/data/diogo/$PASTA" \
        ":neopz-master/Projects/SlopeMonteCarlo/cluster/logs" "$DEST/" || {
        c=$?  # 24: arquivos sumiram durante a copia; 23: alguns arquivos nao foram lidos (ficam as versoes do PC)
        [ $c -eq 23 ] || [ $c -eq 24 ] || { echo "rsync falhou (codigo $c)" >&2; exit $c; }
        echo "AVISO: rsync com codigo $c (veja as mensagens acima); a tabela usa o que ha no PC" >&2
    }
    # partes vazias: o /mnt/data do fermi perdeu o conteudo (2026-11-12); a copia de backup_partes as substitui
    for f in "$DEST/$PASTA"/*/*/partes/*.csv; do
        [ -s "$f" ] || [ ! -e "$f" ] && continue
        b="$DEST/backup_partes/${f#"$DEST/$PASTA"/}"
        if [ -s "$b" ]; then cp "$b" "$f"; echo "AVISO: $f veio vazio do fermi; usei $b" >&2
        else echo "AVISO: $f veio vazio do fermi e nao ha backup" >&2; fi
    done
fi

echo "== tabela (junta as partes e compara com o artigo) =="
"$BIN" resumo "$DEST/$PASTA"

REF="$DEST/$PASTA/ceron_mc_fase1/f50_I0_iii/mc_resumo.csv"
N=$(( $(wc -l < "$REF") - 1 ))
if [ -n "${SEM_BAIXAR:-}" ] || [ -n "$FILA" ]; then SAIDA="$DEST/relatorio_parcial$SUF"; ROT="$N realizacoes (parcial)"
else SAIDA="$DEST/relatorio$SUF"; ROT="$N realizacoes"; fi

echo "== relatorio ($N realizacoes na referencia) =="
(cd "$SCRIPTS" && python3 -I relatorio_sat.py "$DEST/$PASTA" ceron_mc - "$SAIDA" "g15, $NREF niveis h, $ROT") \
    > "$SAIDA.log" 2>&1 || { tail -n 20 "$SAIDA.log"; exit 1; }
echo "pronto: $SAIDA/relatorio.pdf  (figuras em $SAIDA/figs, log em $SAIDA.log)"
