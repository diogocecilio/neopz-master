#!/bin/bash
# status_mc.sh [projeto] [--resumo]: o estado do Monte Carlo submetido por submeter_mc.sh (le
# $SAIDA/<projeto>/campanha.txt; sem projeto, todos os que tem campanha em $SAIDA): a compilacao no no de login, os jobs
# verif-mc-<V2|V4> no SLURM (estado, tempo, no ou motivo da espera) e, por entrada de CASOS, o FS minimo do caso
# deterministico, as realizacoes concluidas, rodando, esperando (as interrompidas, com checkpoint, entre elas), que
# falharam e o total; o tempo medio por realizacao, o FS minimo medio, P_f = P(FS_min < 1) das concluidas ate agora e
# uma estimativa do fim; os avisos, e o comando rsync para baixar os resultados. --resumo mostra tambem resumo_mc.sh de
# cada entrada. Uma realizacao concluida com situacao em SITUACOES_RUPTURA (padrao: sem_equilibrio_geostatico e
# fs_abaixo_de_F0, sem equilibrio nem com a resistencia minima da busca) conta como ruptura; com outra situacao != ok
# (o calculo nao terminou: permanente nao atingido, ...) fica fora de P_f e e contada a parte.
set -u
MC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$MC/.." && pwd)"
cd "$PKG" || exit 1
export LC_ALL=C
source "$PKG/comum.sh"
source "$MC/projetos.sh"
RESUMO=0 PEDIDO="${PROJETO:-}"
for a in "$@"; do
    case "$a" in --resumo) RESUMO=1 ;; *) PEDIDO="$a" ;; esac
done
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
SAIDA="${SAIDA:-$PKG/resultados_mc}"
EU="$VERIF_EU" HOST="$VERIF_MAQ" BOOT="$VERIF_BOOT"
RUP=",${SITUACOES_RUPTURA:-sem_equilibrio_geostatico,fs_abaixo_de_F0},"
# compilacao no no de login (antes de tudo: na primeira vez ainda nao ha jobs)
estado_compilacao
if [ -n "$PEDIDO" ]; then
    mc_projeto "$PEDIDO" > /dev/null || exit 2
    PROJS=("$MC_PROJETO")
else
    PROJS=()
    for p in $(mc_projetos); do [ -f "$SAIDA/$p/campanha.txt" ] && PROJS+=("$p"); done
    if [ ${#PROJS[@]} -eq 0 ]; then
        echo "nenhuma campanha submetida ainda em $SAIDA (rode mc/submeter_mc.sh <projeto> <nos>, ou defina SAIDA)"
        exit 1
    fi
fi

# jobs: "<id> <estado>" de todos os jobs do usuario (para saber se o dono de cada realizacao esta vivo)
JOBS="$(squeue -h -u "$EU" -o '%i %T' 2>/dev/null)"
SQ_OK=$?
command -v squeue > /dev/null 2>&1 || SQ_OK=1
# estado_dono <r/.ativo>: vivo, acabando ou morto (sem trava: morto)
estado_dono() {
    local a="$1" h p j b t
    [ -d "$a" ] || { echo morto; return; }
    if ! read -r h p j b t 2>/dev/null < "$a/dono"; then echo vivo; return; fi
    if [[ "$j" =~ ^[0-9]+([.][0-9]+)?$ ]] && [ $SQ_OK -eq 0 ]; then
        echo "$JOBS" | awk -v j="${j%%.*}" '$1 "" == j "" { print ($2 == "COMPLETING" ? "acabando" : "vivo"); f = 1 } END { if (!f) print "morto" }'
    elif [[ "$j" =~ ^[0-9]+([.][0-9]+)?$ ]]; then echo vivo
    elif [ "$h" = "$HOST" ]; then
        if [ "$b" = "$BOOT" ] && kill -0 "$p" 2>/dev/null; then
            case "$(sed -n 's/^.*) \([A-Za-z]\).*/\1/p' "/proc/$p/stat" 2>/dev/null)" in Z | X | x) echo morto ;; *) echo vivo ;; esac
        else echo morto
        fi
    elif [ -n "$(find "$(dirname "$a")" -mindepth 1 -maxdepth 1 ! -name .ativo ! -name .tomada -mmin -360 | head -n 1)" ]; then echo vivo
    else echo morto
    fi
}
# campo <chave> <arquivo>: o valor de "chave=valor" num resultado.txt
campo() { tr ' ' '\n' < "$2" 2>/dev/null | sed -n "s/^$1=//p" | head -n 1; }

for PJ in "${PROJS[@]}"; do
    mc_projeto "$PJ" || continue
    D="$SAIDA/$MC_PROJETO"
    echo "== $MC_PROJETO"
    E="$PKG/bin/$MC_EXE"
    [ -x "$E" ] && echo "bin/$MC_EXE: $(date -r "$E" '+%F %T') ($(versao_de "$E"))" || echo "bin/$MC_EXE: ainda nao compilado"
    if ! [ -f "$D/campanha.txt" ]; then
        echo "nenhuma campanha submetida ainda em $D (rode mc/submeter_mc.sh $MC_ID <nos>, ou defina SAIDA)"
        continue
    fi
    CASOS="$(sed -n 's/^CASOS="\(.*\)"$/\1/p' "$D/campanha.txt")"
    PRIMEIRA="$(sed -n 's/^PRIMEIRA=//p' "$D/campanha.txt")" ULTIMA="$(sed -n 's/^ULTIMA=//p' "$D/campanha.txt")"
    DET="$(sed -n 's/^DET=//p' "$D/campanha.txt")"
    PRIMEIRA="${PRIMEIRA:-$MC_PRIMEIRA}" ULTIMA="${ULTIMA:-$MC_ULTIMA}" DET="${DET:-1}"
    if [ $SQ_OK -eq 0 ]; then
        n="$(squeue -h -u "$EU" -n "verif-mc-$MC_ID" -o %i 2>/dev/null | wc -l)"
        echo "jobs verif-mc-$MC_ID: $n"
        squeue -h -u "$EU" -n "verif-mc-$MC_ID" -o '%i|%T|%M|%L|%N|%r' 2>/dev/null | awk -F'|' '{
            if ($2 == "RUNNING") printf "  %s RUNNING ha %s (limite em %s) em %s\n", $1, $3, $4, $5
            else printf "  %s %s (%s)\n", $1, $2, $6 }'
    else
        echo "squeue nao respondeu (ou nao existe): o estado dos donos das realizacoes vem so dos arquivos"
    fi
    det_c=0 det_f=0 det_t=0 tot_inc=0 tot_c=0 tot_r=0 tot_e=0 tot_f=0 tot_t=0 soma_min=0 nmin=0
    printf '%-26s %8s %7s %7s %7s %9s %6s %6s %8s %8s %7s\n' entrada "FS det" concl rodando espera "(interr.)" falhou \
        total "min/real" "FS_min" P_f
    while IFS=$'\t' read -r sc DD0 DD; do
        # FS minimo do deterministico
        fsdet="-"
        if [ "$DD" != - ]; then
            rd="$DD/r00000"
            if [ -s "$rd/resultado.txt" ]; then
                fsdet="$(campo FS_min "$rd/resultado.txt")"
                [ "$(campo situacao "$rd/resultado.txt")" = ok ] || fsdet="${fsdet}*"
                det_c=$((det_c + 1))
            elif [ -e "$rd/falhou.txt" ]; then fsdet=falhou; det_f=$((det_f + 1))
            elif [ "$(estado_dono "$rd/.ativo")" = vivo ]; then fsdet=rodando
            else fsdet=espera
            fi
            det_t=$((det_t + 1))
        fi
        # donos: "<realizacao> <estado>" (vivo, acabando ou morto)
        donos="$(for a in "$DD0"/r[0-9]*/.ativo; do
            [ -d "$a" ] || continue
            echo "$(basename "$(dirname "$a")") $(estado_dono "$a")"
        done)"
        linha="$( {
            echo "$donos" | awk 'NF == 2 { print "dono", $1, $2 }'
            find "$DD0" -mindepth 2 -maxdepth 2 \( \( -name resultado.txt -size +0 \) -o -name falhou.txt -o -name checkpoint.bin \) \
                -printf 'arq %P\n' 2>/dev/null
            find "$DD0" -mindepth 2 -maxdepth 2 -name resultado.txt -size +0 -exec cat {} + 2>/dev/null | awk '{ print "res", $0 }'
        } | awk -v a="$PRIMEIRA" -v b="$ULTIMA" -v rup="$RUP" '
            $1 == "dono" { n = substr($2, 2) + 0; if (n >= a && n <= b) dono[n] = $3; next }
            $1 == "arq" { split($2, p, "/"); n = substr(p[1], 2) + 0; if (n < a || n > b) next
                          if (p[2] == "resultado.txt") ok[n] = 1; else if (p[2] == "falhou.txt") fal[n] = 1; else ck[n] = 1; next }
            $1 == "res" { delete v
                          for (i = 2; i <= NF; i++) { k = $i; sub(/=.*/, "", k); x = $i; sub(/^[^=]*=/, "", x); v[k] = x }
                          n = v["real"] + 0; if (n < a || n > b) next
                          if (v["seg"] != "") { s += v["seg"]; ns++ }
                          sit = v["situacao"]; r = (sit != "" && index(rup, "," sit ",") > 0)
                          if (sit != "" && sit != "ok" && !r) { inc++; next }
                          if (v["FS_min"] != "") { F = v["FS_min"] + 0; if (r && F > 0) F = -1
                              nf++; if (F < 1) pf++; if (F > 0) { sf += F; nfs++ } }
                          next }
            END {
                for (n = a; n <= b; n++) {
                    if (n in ok) c++
                    else if (n in fal) f++
                    else if ((n in dono) && dono[n] != "morto") r++
                    else { e++; if (n in ck) it++ }
                }
                printf "%d %d %d %d %d %d %.1f %d %s %s %d\n", c, r, e, it, f, b - a + 1, ns ? s / ns / 60 : 0, ns,
                       nfs ? sprintf("%.4f", sf / nfs) : "-", nf ? sprintf("%.4f", pf / nf) : "-", inc
            }')"
        read -r c r e it f t mmin nmn fsm pf inc <<< "$linha"
        tot_inc=$((tot_inc + inc))
        printf '%-26s %8s %7d %7d %7d %9s %6d %6d %8s %8s %7s\n' "$(basename "$DD0")" "$fsdet" "$c" "$r" "$e" "($it)" "$f" \
            "$t" "$([ "$nmn" -gt 0 ] && echo "$mmin" || echo -)" "$fsm" "$pf"
        tot_c=$((tot_c + c)) tot_r=$((tot_r + r)) tot_e=$((tot_e + e)) tot_f=$((tot_f + f)) tot_t=$((tot_t + t))
        if [ "$nmn" -gt 0 ]; then
            soma_min="$(awk -v s="$soma_min" -v m="$mmin" -v k="$nmn" 'BEGIN { print s + m * k }')"
            nmin=$((nmin + nmn))
        fi
    done < <(PROJETO="$MC_PROJETO" CASOS="$CASOS" DET="$DET" SAIDA="$SAIDA" "$MC/rodar_mc.sh" --pastas)
    media="$(awk -v s="$soma_min" -v k="$nmin" 'BEGIN { if (k) printf "%.1f", s / k; else print "-" }')"
    echo "total: $tot_c de $tot_t concluidas, $tot_r rodando, $tot_e esperando, $tot_f falharam; $media min por realizacao"
    [ "$det_t" -gt 0 ] && echo "deterministicos: $det_c de $det_t concluidos, $det_f falharam (* no FS det: situacao != ok)"
    [ "$tot_inc" -gt 0 ] && echo "$tot_inc realizacoes concluidas com o calculo incompleto (situacao != ok e fora de" \
        "SITUACOES_RUPTURA): fora do FS_min e de P_f (veja o run.err delas)"
    ntent="$(find "$D" -mindepth 3 -maxdepth 3 -name falhas 2>/dev/null | wc -l)"
    [ "$ntent" -gt 0 ] && echo "$ntent realizacoes ja terminaram com erro pelo menos uma vez (arquivo falhas; as com" \
        "falhou.txt nao sao mais tentadas: para tentar de novo, find $D -name falhou.txt -delete; find $D -name falhas -delete)"
    if [ "$tot_r" -gt 0 ] && [ "$nmin" -gt 0 ]; then
        awk -v f="$((tot_t - tot_c - tot_f))" -v m="$media" -v r="$tot_r" 'BEGIN {
            h = f * m / r / 60; printf "fim estimado em %.1f h (%.1f dias), com %d realizacoes ao mesmo tempo\n", h, h / 24, r }'
    fi
    [ -e "$SAIDA/PARAR" ] && echo "ATENCAO: $SAIDA/PARAR existe (todos os projetos): nenhuma realizacao nova comeca;" \
        "para retomar: rm $SAIDA/PARAR e mc/submeter_mc.sh $MC_ID <nos>"
    [ -e "$D/PARAR" ] && echo "ATENCAO: $D/PARAR existe: nenhuma realizacao nova comeca, e os jobs terminam; para" \
        "retomar: mc/submeter_mc.sh $MC_ID <nos>"
    if [ -s "$D/AVISOS" ]; then
        echo "avisos ($D/AVISOS, os ultimos; apague o arquivo depois de ler):"
        tail -n 8 "$D/AVISOS" | cut -c1-250 | sed 's/^/  /'
    fi
    if [ -s "$D/CONTINUACAO_PENDENTE" ]; then
        echo "ATENCAO: o SLURM recusou a continuacao de um job; quando houver vaga na cota, rode mc/submeter_mc.sh $MC_ID <nos>:"
        sed 's/^/  /' "$D/CONTINUACAO_PENDENTE"
    fi
    echo "resultados em $D; para baixar, rode NO PC (pode ser no meio da rodada, so vem o que mudou):"
    echo "  rsync -av --exclude checkpoint.bin --exclude .ativo --exclude .tomada --exclude '.parar.*' --exclude /PARAR" \
        "--exclude /execucoes.txt --exclude /campanha.txt $EU@fermi.cesup.ufrgs.br:$(printf %q "$(cd "$D" && pwd -P)/")" \
        "~/projects/runs/monte_carlo/$MC_PROJETO/"
    if [ "$RESUMO" = 1 ]; then
        PROJETO="$MC_PROJETO" CASOS="$CASOS" DET=0 SAIDA="$SAIDA" "$MC/rodar_mc.sh" --pastas | while IFS=$'\t' read -r sc DD0 DD; do
            [ -d "$DD0" ] && "$MC/resumo_mc.sh" "$DD0"
        done
    fi
done
