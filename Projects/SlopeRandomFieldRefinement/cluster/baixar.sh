#!/bin/bash
# baixar.sh (roda no PC, com a VPN da UFRGS ligada): mostra os jobs da campanha no fermi, baixa a pasta dela em $DATA e
# os logs do SLURM, mostra a tabela por caso e faz o pos-processamento. Pede a senha do fermi uma vez so.
#   ./baixar.sh [destino]                 chuva: destino padrao ~/projects/runs/fermi_chuva; relatorio PDF (posproc/)
#   CAMPANHA=ceron ./baixar.sh [destino]  Ceron: ~/projects/runs/fermi_ceron_ai; P_f de cada amostragem por
#                                         importancia contra o Monte Carlo direto (<destino>/resultado.txt)
#   CAMPANHA=chuva_ai ./baixar.sh         o mesmo para a chuva com amostragem por importancia (~/projects/runs/fermi_chuva_ai)
#   CAMPANHA=etapa1 ./baixar.sh           etapas 0 e 1 (~/projects/runs/fermi_etapa1): tambem as comparacoes pareadas (par=)
#                                         e os FORM
#   SEM_BAIXAR=1 ./baixar.sh     so o pos-processamento, com os dados ja baixados (sem VPN)
# Chuva: o relatorio vai para <destino>/relatorio_parcial enquanto ha jobs na fila e para <destino>/relatorio depois.
# Ambiente: CAMPANHA, FERMI (padrao diogo@fermi.cesup.ufrgs.br), REMOTO (padrao /mnt/data/diogo/<pasta da campanha>).
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$AQUI/campanha.sh"
DEST="${1:-$HOME/projects/runs/$DESTINO}"
HOST="${FERMI:-diogo@fermi.cesup.ufrgs.br}"
REMOTO="${REMOTO:-/mnt/data/diogo/$PASTA}"
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
    FILA="$("${SSH[@]}" "$HOST" "squeue -h -u \$USER -n $NOME -o '%.8i %.2t %.12M %R'")"
    if [ -n "$FILA" ]; then echo "$FILA"; else echo "nenhum job $NOME na fila"; fi
    echo "== tabela no fermi =="
    "${SSH[@]}" "$HOST" "CAMPANHA=$CAMPANHA SAIDA=$REMOTO neopz-master/Projects/SlopeRandomFieldRefinement/cluster/status.sh" |
        tail -n +2 || true

    echo "== baixando para $DEST =="
    rsync -az --info=progress2 -e "${SSH[*]}" --exclude .tarefas/ "$HOST:$REMOTO" \
        ":neopz-master/Projects/SlopeRandomFieldRefinement/cluster/logs" "$DEST/" || {
        c=$? # 24: arquivos sumiram durante a copia; 23: alguns arquivos nao foram lidos (ficam as versoes do PC)
        [ $c -eq 23 ] || [ $c -eq 24 ] || { echo "rsync falhou (codigo $c)" >&2; exit $c; }
        echo "AVISO: rsync com codigo $c (veja as mensagens acima); o relatorio usa o que ha no PC" >&2
    }
fi

if [ "$CAMPANHA" != chuva ]; then
    # o resumo instante a instante de cada caso (Monte Carlo direto); cada amostragem por importancia contra o Monte Carlo
    # direto de referencia (o piloto=, ou o ref= das opcoes da fila); cada caso com par= contra o seu par, realizacao a
    # realizacao (compara_pares.py); e os FORM (form_M*, do form.sh)
    R="$DEST/resultado.txt" D="$DEST/$(basename "$REMOTO")"
    python3 -I -B "$POS/resumo_tempos.py" "$D" "$AQUI/$CASOS_TXT" > "$R" 2>&1 || true
    echo >> "$R"
    while read -r id fase n opc; do
        [[ -z "$id" || "$id" == \#* ]] && continue
        [ -d "$D/$id" ] || continue
        ref="" par=""
        [[ " $opc " =~ \ piloto=@SAIDA@/([^ ]+)\  ]] && ref="${BASH_REMATCH[1]}"
        [[ " $opc " =~ \ ref=([^ ]+)\  ]] && ref="${BASH_REMATCH[1]}"
        [[ " $opc " =~ \ par=([^ ]+)\  ]] && par="${BASH_REMATCH[1]}"
        if [ -n "$ref" ] && [ -d "$D/$ref" ]; then
            { echo "== $id (Monte Carlo direto de referencia: $ref)"
              python3 -I -B "$POS/amostragem_importancia.py" estima "$D/$id" "$D/$ref" 2>&1 || echo "(sem estimativa ainda)"
              echo; } >> "$R"
        fi
        if [ -n "$par" ] && [ -d "$D/$par" ]; then
            { echo "== $id contra $par, realizacao a realizacao"
              python3 -I -B "$POS/compara_pares.py" "$D/$par" "$D/$id" 2>&1 || echo "(sem comparacao ainda)"
              echo; } >> "$R"
        fi
    done < "$AQUI/$CASOS_TXT"
    if ls -d "$D"/form_M* > /dev/null 2>&1; then
        { echo "== FORM"; python3 -I -B "$POS/resumo_form.py" "$D"/form_M* 2>&1; echo; } >> "$R"
    fi
    cat "$R"
    echo "pronto: $R"
    exit 0
fi
if [ -n "${SEM_BAIXAR:-}" ] || [ -n "$FILA" ]; then SAIDA="$DEST/relatorio_parcial"; else SAIDA="$DEST/relatorio"; fi
mkdir -p "$SAIDA"
echo "== estatisticas e relatorio em $SAIDA =="
python3 -I -B "$POS/dados_chuva.py" "$DEST/$(basename "$REMOTO")" "$SAIDA/estatisticas.json" > "$SAIDA.log" 2>&1 &&
    python3 -I -B "$POS/relatorio_chuva.py" "$SAIDA/estatisticas.json" "$SAIDA" >> "$SAIDA.log" 2>&1 ||
    { tail -n 20 "$SAIDA.log"; exit 1; }
grep -v "^aviso" "$SAIDA.log" | head -n 40
echo "pronto: $SAIDA/relatorio.pdf  (figuras em $SAIDA/figs, log em $SAIDA.log)"
