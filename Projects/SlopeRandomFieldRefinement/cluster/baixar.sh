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
    # realizacao (compara_pares.py); e os FORM (form_*, do form.sh)
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
    if ls -d "$D"/form_* > /dev/null 2>&1; then
        { echo "== FORM"; python3 -I -B "$POS/resumo_form.py" "$D"/form_* 2>&1; echo; } >> "$R"
    fi
    # etapa 3b: quantas realizacoes a ISMC precisa para o CoV alvo, por centro e por r, contra os diretos das etapas 1 e 3
    if [ "$CAMPANHA" = etapa3b ]; then
        { echo "== ISMC com parada pelo CoV (Projects/Cho2010/ismc_parada.py)"
          python3 -I -B "$AQUI/../../Cho2010/ismc_parada.py" "$D" "$AQUI/$CASOS_TXT" 2>&1; echo; } >> "$R"
    fi
    # etapa3b_gim: cada caso GIM contra a sua referencia SRM, realizacao a realizacao (Projects/Cho2010/pareado_gim.py;
    # as referencias sao as ja baixadas das etapas 1, 3 e 3b), e o line sampling GIM contra o SRM, linha a linha
    if [ "$CAMPANHA" = etapa3b_gim ]; then
        RUNS="$HOME/projects/runs" PAR="$AQUI/../../Cho2010/pareado_gim.py"
        { for c in conf_srm_r05:fermi_etapa1/etapa1/r05_direto conf_srm_r07:fermi_etapa3/etapa3/r07_direto; do
              id=${c%%:*} ref=${c#*:}
              [ -d "$D/$id" ] || continue
              echo "== conferencia de versao: $id contra $ref (FS iguais)"
              python3 -I -B "$POS/compara_pares.py" "$RUNS/$ref" "$D/$id" 2>&1 | head -n 12 || true
              echo
          done
          for c in r05_nys_M20_s21: r05_nys_M100_s21: r05_nys_gim20_s21: r07_nys_M20_s21: r07_nys_M100_s21: r07_nys_gim20_s21:; do
              id=${c%%:*} ref=${c#*:}
              [ -n "$ref" ] || ref="fermi_etapa3b/etapa3b/$id"
              [ -d "$D/$id" ] || continue
              echo "== $id: SRM x GIM pareado"
              python3 -I -B "$PAR" mc "$RUNS/$ref" "$D/$id" 2>&1 || true
              echo
          done
          for c in r05_M20 r07_M20 r05_gim2000 r07_gim2000; do
              [ -f "$D/ls/$c/linhas.csv" ] || continue
              echo "== line sampling $c: SRM x GIM linha a linha"
              python3 -I -B "$PAR" ls "$RUNS/etapa3b_ls/ls/$c" "$D/ls/$c" 2>&1 || true
              echo
          done
          # N alvo com cov_n=500 (a regra do SRM na etapa 3b; a fila rodou com 2000 so para parear), o custo com os FORM
          # da etapa 3b (os centros sao os de la) e o z contra o direto SRM (sem direto GIM: decisao do usuario)
          echo "== ISMC GIM: parada pelo CoV com cov_n=500 (como o SRM), FORM da etapa 3b, z contra o direto SRM"
          python3 -I -B "$AQUI/../../Cho2010/ismc_parada.py" "$D" "$AQUI/$CASOS_TXT" cov_n=500 \
              forms="$RUNS/fermi_etapa3b/etapa3b" r05="$RUNS/fermi_etapa1/etapa1/r05_direto" r07="$RUNS/fermi_etapa3/etapa3/r07_direto" 2>&1 || true
          echo; } >> "$R"
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
