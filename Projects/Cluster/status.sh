#!/bin/bash
# status.sh: o estado de tudo o que este pacote roda no cluster: a compilacao em segundo plano, os executaveis de bin/
# (com a versao), os jobs verif-* no SLURM, os avisos (logs/avisos.txt), os casos completos em $SAIDA_CASOS (padrao
# $DATA/verificacao: por projeto, os concluidos com a situacao e o tempo, os rodando e os sem resultado), o ultimo job
# de testes (logs/testes_<job>/resumo.txt) e uma linha por campanha de Monte Carlo (os detalhes, o fim estimado e o
# rsync: mc/status_mc.sh <projeto>). No fim, o comando rsync para baixar os casos completos para o PC.
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PKG" || exit 1
export LC_ALL=C
source "$PKG/comum.sh"
source "$PKG/mc/projetos.sh"
estado_compilacao
echo "executaveis (bin/):"
for p in $(projetos_verif); do
    e="$PKG/bin/$(alvo_exe "$p")"
    if [ -x "$e" ]; then
        n="$(fonte_nova "$p" "$e")"
        echo "  ${e##*/}: $(versao_de "$e")${n:+ (fontes mais novos: ${n#../../}; recompile)}"
    else
        echo "  ${e##*/}: nao compilado (./compilar.sh)"
    fi
done
if command -v squeue > /dev/null 2>&1; then
    if fila="$(squeue -h -u "$VERIF_EU" -o '%i|%j|%T|%M|%L|%N|%r' 2>&1)"; then
        fila="$(echo "$fila" | awk -F'|' '$2 ~ /^verif-/')"
        echo "jobs verif-*: $(echo -n "$fila" | grep -c .)"
        [ -n "$fila" ] && echo "$fila" | awk -F'|' '{
            if ($3 == "RUNNING") printf "  %s %s RUNNING ha %s (limite em %s) em %s\n", $1, $2, $4, $5, $6
            else printf "  %s %s %s (%s)\n", $1, $2, $3, $7 }'
    else
        echo "squeue nao respondeu: $fila"
    fi
fi
if [ -s "$PKG/logs/avisos.txt" ]; then
    echo "avisos (logs/avisos.txt, os ultimos; apague o arquivo depois de ler):"
    tail -n 6 "$PKG/logs/avisos.txt" | cut -c1-250 | sed 's/^/  /'
fi

# casos completos
SAIDA_CASOS="${SAIDA_CASOS:-${DATA:+$DATA/verificacao}}"
SAIDA_CASOS="${SAIDA_CASOS:-$PKG/resultados_casos}"
if [ -d "$SAIDA_CASOS" ]; then
    echo "casos completos em $SAIDA_CASOS:"
    for p in $(projetos_verif); do
        [ -d "$SAIDA_CASOS/$p" ] || continue
        linhas=""
        for d in "$SAIDA_CASOS/$p"/*/; do
            d="${d%/}"
            c="${d##*/}"
            case "$c" in *.ant_*) continue ;; esac # as rodadas antigas (REFAZER=1)
            [ -f "$d/caso.txt" ] || [ -d "$d/.ativo" ] || continue
            if [ -s "$d/resultado.txt" ]; then
                linhas="$linhas    $c: $(tr ' ' '\n' < "$d/resultado.txt" | grep -E '^(situacao|seg)=' | tr '\n' ' ')\n"
            elif [ -d "$d/.ativo" ]; then
                linhas="$linhas    $c: rodando ($(cat "$d/.ativo/dono" 2>/dev/null | awk '{ print "job " $1 " em " $2 }'))\n"
            else
                linhas="$linhas    $c: SEM resultado.txt ($(tail -n 1 "$d/caso.txt"); veja $d/log.txt)\n"
            fi
        done
        [ -n "$linhas" ] && { echo "  $p:"; printf '%b' "$linhas"; }
    done
fi

# o ultimo job de testes
u="$(ls -dt "$PKG"/logs/testes_*/ 2>/dev/null | head -n 1)"
if [ -n "$u" ] && [ -f "$u/resumo.txt" ]; then
    echo "ultimo teste (${u%/}/resumo.txt):"
    sed 's/^/  /' "$u/resumo.txt"
    grep -q ' fim: ' "$u/resumo.txt" || echo "  (ainda rodando, ou o job acabou antes do fim)"
fi

# Monte Carlo
SAIDA_MC="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
SAIDA_MC="${SAIDA_MC:-$PKG/resultados_mc}"
for p in $(mc_projetos); do
    D="$SAIDA_MC/$p"
    [ -f "$D/campanha.txt" ] || continue
    mc_projeto "$p" > /dev/null || continue
    a="$(sed -n 's/^PRIMEIRA=//p' "$D/campanha.txt")" b="$(sed -n 's/^ULTIMA=//p' "$D/campanha.txt")"
    casos="$(sed -n 's/^CASOS="\(.*\)"$/\1/p' "$D/campanha.txt")"
    # so as pastas e a faixa de realizacoes da campanha atual
    read -r nc nf < <(PROJETO="$p" CASOS="$casos" SAIDA="$SAIDA_MC" "$PKG/mc/rodar_mc.sh" --pastas 2>/dev/null |
        while IFS=$'\t' read -r _ pd _; do
            find "$pd" -mindepth 2 -maxdepth 2 \( \( -name resultado.txt -size +0 \) -o -name falhou.txt \) -printf '%P\n' 2>/dev/null
        done | awk -F/ -v a="$a" -v b="$b" '{ n = substr($1, 2) + 0; if (n < a || n > b) next
            if ($2 == "resultado.txt") c++; else f++ } END { print c + 0, f + 0 }')
    echo "Monte Carlo $p: ${nc:-0} de $(( (b - a + 1) * $(echo $casos | wc -w) )) realizacoes concluidas (${nf:-0} falharam);" \
        "detalhes: mc/status_mc.sh $MC_ID$([ -e "$D/PARAR" ] || [ -e "$SAIDA_MC/PARAR" ] && echo '; PARADO (PARAR)')"
done
if [ -d "$SAIDA_CASOS" ]; then
    echo "para baixar os casos completos, rode NO PC:"
    echo "  rsync -av --exclude .ativo --exclude /testes/ $VERIF_EU@fermi.cesup.ufrgs.br:$(printf %q "$(cd "$SAIDA_CASOS" && pwd -P)/")" \
        "~/projects/runs/verificacao/"
fi
