#!/bin/bash
# junta_mc.sh <pasta da combinacao> [...]: junta os CSV das realizacoes concluidas (pastas r<nnnnn>/ com resultado.txt):
#   mc_resumo.csv     <- resumo.csv de cada realizacao (uma linha por realizacao)
#   mc_fs.csv         <- fs.csv (uma linha por realizacao e instante do FS)
#   mc_passos.csv     <- passos.csv (uma linha por realizacao e passo do transiente), so com TUDO=1
#   mc_tentativas.csv <- fs_tentativas.csv (cada solucao da busca do SRM), so com TUDO=1
# A primeira coluna e a realizacao (o numero da pasta). As colunas sao a uniao das colunas de todas as realizacoes: as de
# uma versao anterior do programa (sem resumo.csv, com fs.csv de 3 colunas) entram com o que tiverem, e o resto vazio.
# Refeito do zero a partir das pastas a cada chamada (troca atomica), entao pode rodar durante o Monte Carlo;
# rodar_mc.sh, resumo_mc.sh e plot_mc.py o chamam.
set -u
# junta <saida> <arquivos...>: dois passes do awk, o 1o so le os cabecalhos (uniao das
# colunas, na ordem em que aparecem) e o 2o escreve cada linha nessa ordem
junta() {
    local saida="$1"
    shift
    [ $# -gt 0 ] || { rm -f "$saida"; return 0; }
    awk -F, '
        passo == 1 {
            if (FNR == 1) { n = split($0, h, ","); for (i = 1; i <= n; i++) if (h[i] != "real" && !(h[i] in pos)) { pos[h[i]] = ++nu; uni[nu] = h[i] } }
            nextfile
        }
        !impresso { printf "real"; for (j = 1; j <= nu; j++) printf ",%s", uni[j]; printf "\n"; impresso = 1 }
        FNR == 1 { d = FILENAME; sub(/\/[^\/]*$/, "", d); r = d; sub(/.*\/r/, "", r); r += 0; nc = split($0, nome, ","); next }
        $0 != "" {
            nv = split($0, v, ",")
            for (j = 1; j <= nu; j++) o[j] = ""
            for (i = 1; i <= nc && i <= nv; i++) if (nome[i] in pos) o[pos[nome[i]]] = v[i]
            printf "%d", r; for (j = 1; j <= nu; j++) printf ",%s", o[j]; printf "\n"
        }' passo=1 "$@" passo=2 "$@" > "$saida.tmp" && mv -f "$saida.tmp" "$saida"
}
for D in "$@"; do
    [ -d "$D" ] || { echo "$D: nao e uma pasta" >&2; continue; }
    tmp="$(mktemp -d)"
    res=() fs=() pas=() ten=()
    for r in "$D"/r[0-9]*/; do
        r="${r%/}"
        [ -s "$r/resultado.txt" ] || continue # em andamento ou interrompida
        if [ -s "$r/resumo.csv" ]; then
            res+=("$r/resumo.csv")
        else # versao anterior: o resumo vem de resultado.txt (chave=valor)
            mkdir -p "$tmp/$(basename "$r")"
            tr ' ' '\n' < "$r/resultado.txt" | awk -F= 'NF == 2 { k[++n] = $1; v[n] = $2 } END {
                for (i = 1; i <= n; i++) printf "%s%s", (i > 1 ? "," : ""), k[i]; printf "\n"
                for (i = 1; i <= n; i++) printf "%s%s", (i > 1 ? "," : ""), v[i]; printf "\n" }' > "$tmp/$(basename "$r")/resumo.csv"
            res+=("$tmp/$(basename "$r")/resumo.csv")
        fi
        [ -s "$r/fs.csv" ] && fs+=("$r/fs.csv")
        [ -s "$r/passos.csv" ] && pas+=("$r/passos.csv")
        [ -s "$r/fs_tentativas.csv" ] && ten+=("$r/fs_tentativas.csv")
    done
    junta "$D/mc_resumo.csv" ${res[@]+"${res[@]}"}
    junta "$D/mc_fs.csv" ${fs[@]+"${fs[@]}"}
    if [ "${TUDO:-0}" = 1 ]; then
        junta "$D/mc_passos.csv" ${pas[@]+"${pas[@]}"}
        junta "$D/mc_tentativas.csv" ${ten[@]+"${ten[@]}"}
    fi
    rm -rf "$tmp"
    echo "$D: ${#res[@]} realizacoes em mc_resumo.csv e mc_fs.csv$([ "${TUDO:-0}" = 1 ] && echo ', mc_passos.csv e mc_tentativas.csv')"
done
