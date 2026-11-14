#!/bin/bash
# status.sh: os jobs da campanha na fila, as ultimas linhas dos logs dos processos e, por caso, as tarefas feitas, rodando
# e com erro, as realizacoes ja calculadas e o P_f ate agora (so com bash e awk: roda no fermi). Cada realizacao conta com
# o menor FS_min dos seus instantes (na chuva com 0 e 24 h, o de 24 h; no Ceron, o unico) e so se todos tiverem situacao
# ok; falha: FS < 1. Com a coluna log_peso (amostragem por importancia), P_f = soma(peso * falha) / N, e o CoV vem da
# amostra. CAMPANHA e SAIDA como em submeter.sh.
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./campanha.sh
SAIDA="${SAIDA:-${DATA:+$DATA/$PASTA}}"
SAIDA="${SAIDA:-$PWD/resultados_$CAMPANHA}"
command -v squeue > /dev/null && squeue -u "${USER:-$(id -un)}" -n "$NOME"
[ -f "$SAIDA/campanha.txt" ] || { echo "nada em $SAIDA ainda"; exit 0; }
read -r IREF BLOCO < "$SAIDA/campanha.txt"
echo "campanha $CAMPANHA em $SAIDA: IREF = $IREF, blocos de $BLOCO"
if [ -d "$SAIDA/logs" ]; then
    echo "ultimas linhas dos logs:"
    for f in "$SAIDA"/logs/p*.txt; do [ -e "$f" ] && echo "  $(basename "$f" .txt): $(tail -n 1 "$f")"; done | sort -V | tail -n 6
    echo "erros: $(cat "$SAIDA"/logs/p*.txt 2> /dev/null | grep -c ERRO)"
fi
printf "\n%-14s %7s %7s %5s %7s %7s %8s %7s\n" caso tarefas rodando erros real falhas "P_f(%)" "CoV(%)"
while read -r id fase n opc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    d="$SAIDA/$id"
    [ -d "$d" ] || continue
    feitas=$(ls "$d"/mc_mapa_*.csv 2> /dev/null | wc -l)
    travas=$(ls "$d"/.tarefas/ 2> /dev/null | wc -l)
    erros=$(grep -ls '^ERRO' "$d"/.tarefas/* 2> /dev/null | wc -l)
    rodando=$((travas - feitas - erros))
    ls "$d"/mc_[0-9]*_[0-9]*.csv > /dev/null 2>&1 || { printf "%-14s %7d %7d %5d\n" "$id" "$feitas" "$rodando" "$erros"; continue; }
    awk -F, -v id="$id" -v feitas="$feitas" -v rodando="$rodando" -v erros="$erros" '
        FNR == 1 { delete c; for (i = 1; i <= NF; i++) c[$i] = i; next }
        {
            r = $c["real"] + 0; x = $c["FS_min"] + 0
            if (r < 0) next
            if (!(r in fs) || x < fs[r]) fs[r] = x
            if ($c["situacao"] != "ok" || x <= 0) ruim[r] = 1
            lw[r] = ("log_peso" in c) ? $c["log_peso"] + 0 : 0
        }
        END {
            for (r in fs) if (!(r in ruim)) { n++; q = (fs[r] < 1) ? exp(lw[r]) : 0; f += (fs[r] < 1); s += q; s2 += q * q }
            p = n ? s / n : 0
            cv = (n > 1 && p > 0) ? 100 * sqrt((s2 - n * p * p) / (n - 1) / n) / p : 0
            printf "%-14s %7d %7d %5d %7d %7d %8.2f %7.1f\n", id, feitas, rodando, erros, n, f, 100 * p, cv
        }' "$d"/mc_[0-9]*_[0-9]*.csv
done < "$CASOS_TXT"
