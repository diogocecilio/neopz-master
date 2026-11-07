#!/bin/bash
# status_mc.sh [--resumo]: o estado do Monte Carlo submetido por submeter_mc.sh (le $SAIDA/campanha.txt): a compilacao
# no no de login, os jobs slopeRF no SLURM (estado, tempo, no ou motivo da espera) e, por combinacao, o FS minimo do
# caso deterministico, as realizacoes concluidas, rodando, esperando (as interrompidas, com checkpoint, entre elas),
# que falharam e o total; o tempo medio por realizacao, o FS minimo medio, P_f = P(FS_min < 1) das concluidas ate agora e
# uma estimativa do fim. --resumo mostra tambem resumo_mc.sh de cada combinacao (P_f por instante).
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PKG" || exit 1
export LC_ALL=C
EU="${USER:-$(id -un)}"
HOST="$(hostname -s 2>/dev/null || uname -n)"
HOST="${HOST%%.*}"
HOST="${HOST:-sem_nome}"
BOOT="$(cat /proc/sys/kernel/random/boot_id 2>/dev/null || echo -)"
# compilacao no no de login (antes de tudo: na primeira vez ainda nao ha jobs)
if [ -f logs/compilar_mc.pid ]; then
    read -r cpid cno cboot cini < logs/compilar_mc.pid
    if [ "$cno" = "$HOST" ] && [ "$cboot" = "$BOOT" ] && kill -0 "$cpid" 2>/dev/null &&
        tr '\0' ' ' 2>/dev/null < "/proc/$cpid/cmdline" | grep -q compilar.sh; then
        echo "compilacao: rodando desde $cini (pid $cpid): $(grep -o '\[ *[0-9]*%\]' logs/compilar_mc.log | tail -n 1); ao fim submete: $(cat logs/compilar_mc.pedido 2>/dev/null)"
    elif [ -f logs/compilar_mc.log ] && ! grep -q "^Compilado:" logs/compilar_mc.log; then
        echo "compilacao: terminou SEM gerar o executavel (comecou em $cini); veja: tail -30 logs/compilar_mc.log"
    elif [ -f logs/compilar_mc.log ] && grep -q "recusou\|nenhum job submetido\|nao submetido" logs/compilar_mc.log; then
        echo "compilacao: terminou, mas a submissao dos jobs falhou; veja: tail -5 logs/compilar_mc.log"
    fi
fi
[ -x bin/SlopeRandomField ] && echo "bin/SlopeRandomField: $(date -r bin/SlopeRandomField '+%F %T') ($(cat bin/versao.txt 2>/dev/null || echo 'versao nao registrada'))"
SAIDA="${SAIDA:-${DATA:+$DATA/monte_carlo}}"
SAIDA="${SAIDA:-$PKG/resultados_mc}"
[ -f "$SAIDA/campanha.txt" ] || { echo "nenhuma campanha submetida ainda em $SAIDA (rode ./submeter_mc.sh <nos>, ou defina SAIDA)"; exit 1; }
COMBS="" PRIMEIRA=0 ULTIMA=999 DET=1 VARIANTE=iii
# shellcheck disable=SC1091
source "$SAIDA/campanha.txt"
COMBS="$(echo $COMBS | sed -E 's/,(f10|f50|f100):/ \1:/g')"

# jobs: "<id> <estado>" de todos os jobs do usuario (para saber se o dono de cada realizacao esta vivo)
JOBS="$(squeue -h -u "$EU" -o '%i %T' 2>/dev/null)"
SQ_OK=$?
if [ $SQ_OK -eq 0 ]; then
    n="$(squeue -h -u "$EU" -n slopeRF -o %i | wc -l)"
    echo "jobs slopeRF: $n"
    squeue -h -u "$EU" -n slopeRF -o '%i|%T|%M|%L|%N|%r' 2>/dev/null | awk -F'|' '{
        if ($2 == "RUNNING") printf "  %s RUNNING ha %s (limite em %s) em %s\n", $1, $3, $4, $5
        else printf "  %s %s (%s)\n", $1, $2, $6 }'
else
    echo "squeue nao respondeu: o estado dos donos das realizacoes vem so dos arquivos"
fi

# uma linha por combinacao: le, de uma vez, os marcadores de cada realizacao (find) e os resultado.txt
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
det_c=0 det_f=0 det_t=0 tot_inc=0 limite=()
printf '%-24s %8s %7s %7s %7s %9s %6s %6s %8s %8s %7s\n' combinacao "FS det" concl rodando espera "(interr.)" falhou total \
    "min/real" "FS_min" P_f
tot_c=0 tot_r=0 tot_e=0 tot_f=0 tot_t=0 soma_min=0 nmin=0
while read -r sc D DD; do
    # FS minimo do deterministico
    fsdet="-"
    if [ "$DD" != - ]; then
        rd="$DD/r00000"
        if [ -s "$rd/resultado.txt" ]; then
            fsdet="$(sed -n 's/.*FS_min=\([^ ]*\).*/\1/p' "$rd/resultado.txt")"
            grep -q "situacao=ok" "$rd/resultado.txt" || fsdet="${fsdet}*"
            det_c=$((det_c + 1))
        elif [ -e "$rd/falhou.txt" ]; then fsdet=falhou; det_f=$((det_f + 1))
        elif [ "$(estado_dono "$rd/.ativo")" = vivo ]; then fsdet=rodando
        else fsdet=espera
        fi
        det_t=$((det_t + 1))
    fi
    # donos: "<realizacao> <estado>" (vivo, acabando ou morto)
    donos="$(for a in "$D"/r[0-9]*/.ativo; do
        [ -d "$a" ] || continue
        echo "$(basename "$(dirname "$a")") $(estado_dono "$a")"
    done)"
    linha="$( {
        echo "$donos" | awk 'NF == 2 { print "dono", $1, $2 }'
        find "$D" -mindepth 2 -maxdepth 2 \( \( -name resultado.txt -size +0 \) -o -name falhou.txt -o -name checkpoint.bin \) \
            -printf 'arq %P\n' 2>/dev/null
        find "$D" -mindepth 2 -maxdepth 2 -name resultado.txt -size +0 -exec cat {} + 2>/dev/null | awk '{ print "res", $0 }'
    } | awk -v a="$PRIMEIRA" -v b="$ULTIMA" '
        $1 == "dono" { n = substr($2, 2) + 0; if (n >= a && n <= b) dono[n] = $3; next }
        $1 == "arq" { split($2, p, "/"); n = substr(p[1], 2) + 0; if (n < a || n > b) next
                      if (p[2] == "resultado.txt") ok[n] = 1; else if (p[2] == "falhou.txt") fal[n] = 1; else ck[n] = 1; next }
        $1 == "res" { delete v; for (i = 2; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
                      n = v["real"] + 0; if (n < a || n > b) next
                      if (v["seg"] != "") { s += v["seg"]; ns++ }
                      if (v["situacao"] != "" && v["situacao"] != "ok") { inc++; next }
                      if (v["FS_min"] != "") { F = v["FS_min"] + 0; nf++; if (F < 1) pf++; if (F < 1.01) pl++; if (F > 0) { sf += F; nfs++ } } next }
        END {
            for (n = a; n <= b; n++) {
                if (n in ok) c++
                else if (n in fal) f++
                else if ((n in dono) && dono[n] != "morto") r++
                else { e++; if (n in ck) it++ }
            }
            printf "%d %d %d %d %d %d %.1f %d %s %s %d %s\n", c, r, e, it, f, b - a + 1, ns ? s / ns / 60 : 0, ns,
                   nfs ? sprintf("%.4f", sf / nfs) : "-", nf ? sprintf("%.4f", pf / nf) : "-", inc, nf ? sprintf("%.4f", pl / nf) : "-"
        }')"
    read -r c r e it f t mmin nmn fsm pf inc plim <<< "$linha"
    # variantes acopladas com a resistencia real no transiente (F_TRANS = 1, o padrao): o FS do SRM nao cai abaixo de 1
    ft="$(sed -n 's/^F_TRANS=//p' "$D/parametros.txt" 2>/dev/null)"
    # (F_TRANS >= 1: o FS fica >= F_TRANS); as interrompidas (situacao != ok) podem ser rupturas no transiente
    case "$VARIANTE" in i | ii)
        if awk -v v="$ft" 'BEGIN { exit !(v == "" || v + 0 >= 1) }' && [ "$plim" != - ]; then
            limite+=("$(basename "$D") $plim$([ "${inc:-0}" -gt 0 ] && echo " (mais $inc com o transiente interrompido)")")
        fi ;;
    esac
    tot_inc=$((tot_inc + inc))
    printf '%-24s %8s %7d %7d %7d %9s %6d %6d %8s %8s %7s\n' "$(basename "$D")" "$fsdet" "$c" "$r" "$e" "($it)" "$f" "$t" \
        "$([ "$nmn" -gt 0 ] && echo "$mmin" || echo -)" "$fsm" "$pf"
    tot_c=$((tot_c + c)) tot_r=$((tot_r + r)) tot_e=$((tot_e + e)) tot_f=$((tot_f + f)) tot_t=$((tot_t + t))
    [ "$nmn" -gt 0 ] && { soma_min="$(awk -v s="$soma_min" -v m="$mmin" -v k="$nmn" 'BEGIN { print s + m * k }')"; nmin=$((nmin + nmn)); }
done < <(COMBS="$COMBS" DET="$DET" VARIANTE="$VARIANTE" SAIDA="$SAIDA" ./rodar_mc.sh --pastas)
media="$(awk -v s="$soma_min" -v k="$nmin" 'BEGIN { if (k) printf "%.1f", s / k; else print "-" }')"
echo "total: $tot_c de $tot_t concluidas, $tot_r rodando, $tot_e esperando, $tot_f falharam; $media min por realizacao"
[ "$det_t" -gt 0 ] && echo "deterministicos: $det_c de $det_t concluidos, $det_f falharam"
if [ ${#limite[@]} -gt 0 ]; then
    echo "variante $VARIANTE com F_TRANS >= 1: o FS do SRM nao cai abaixo de 1 (o transiente segura o talude no limite), entao"
    echo "  P_f = P(FS_min < 1) subestima a ruptura. Chegaram ao limite, P(FS_min < 1.01): $(printf '%s, ' "${limite[@]}" |
        sed 's/, $//'). Ver F_TRANS no README.md"
fi
[ "$tot_inc" -gt 0 ] && echo "$tot_inc realizacoes concluidas com o transiente interrompido (situacao != ok; * no FS det):" \
    "fora do FS_min e de P_f (veja run.err delas)"
ntent="$(find "$SAIDA" -mindepth 3 -maxdepth 3 -name falhas 2>/dev/null | wc -l)"
[ "$ntent" -gt 0 ] && echo "$ntent realizacoes ja terminaram com erro pelo menos uma vez (arquivo falhas; as com falhou.txt nao sao" \
    "mais tentadas: para tentar de novo, find $SAIDA -name falhou.txt -delete; find $SAIDA -name falhas -delete)"
if [ "$tot_r" -gt 0 ] && [ "$nmin" -gt 0 ]; then
    awk -v f="$((tot_t - tot_c - tot_f))" -v m="$media" -v r="$tot_r" 'BEGIN {
        h = f * m / r / 60; printf "fim estimado em %.0f h (%.1f dias), com %d realizacoes ao mesmo tempo\n", h, h / 24, r }'
fi
[ -e "$SAIDA/PARAR" ] && echo "ATENCAO: $SAIDA/PARAR existe: nenhuma realizacao nova comeca, e os jobs terminam; para retomar: ./submeter_mc.sh <nos>"
if [ -s "$SAIDA/AVISOS" ]; then
    echo "avisos ($SAIDA/AVISOS, os ultimos; apague o arquivo depois de ler):"
    tail -n 8 "$SAIDA/AVISOS" | cut -c1-250 | sed 's/^/  /'
fi
if [ -s "$SAIDA/CONTINUACAO_PENDENTE" ]; then
    echo "ATENCAO: o SLURM recusou a continuacao de um job; quando houver vaga na cota, rode ./submeter_mc.sh <nos>:"
    sed 's/^/  /' "$SAIDA/CONTINUACAO_PENDENTE"
fi
# destino no PC com o nome da pasta de SAIDA: com um destino fixo, uma SAIDA propria (como a da sensibilidade) cairia na
# pasta da campanha principal, e uma pasta de mesmo nome (f50_I36_i) sobrescreveria as realizacoes dela
saida_real="$(cd "$SAIDA" && pwd -P)"
echo "resultados em $SAIDA; para baixar, rode NO PC (pode ser no meio da rodada, so vem o que mudou):"
echo "  rsync -av --exclude checkpoint.bin --exclude .ativo --exclude .tomada --exclude '.parar.*' --exclude /PARAR" \
    "--exclude /execucoes.txt --exclude /campanha.txt $EU@fermi.cesup.ufrgs.br:$(printf %q "$saida_real/")" \
    "~/projects/runs/$(printf %q "$(basename "$saida_real")")/"
if [ "${1:-}" = --resumo ]; then
    COMBS="$COMBS" DET=0 VARIANTE="$VARIANTE" SAIDA="$SAIDA" ./rodar_mc.sh --pastas | while read -r sc D DD; do
        [ -d "$D" ] && ./resumo_mc.sh "$D"
    done
fi
