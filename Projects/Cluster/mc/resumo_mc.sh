#!/bin/bash
# resumo_mc.sh <pasta> [...]: junta as realizacoes concluidas de uma pasta do Monte Carlo (junta_mc.sh) e mostra a
# estatistica: o FS minimo de cada realizacao (media, desvio, P_f = P(FS_min < 1) e o CoV da estimativa,
# sqrt((1 - P_f) / (N P_f)), Vargas Ceron et al. 2025, eq. 68) e a correlacao (Pearson) do FS minimo com as medias dos
# campos que o programa grava (no solo, na banda de ruptura, na massa que desliza); com mais de um instante do FS por
# realizacao (um transiente), antes o mesmo em cada instante. FS <= 0 (sem equilibrio nem com a resistencia minima da
# busca) conta como ruptura e fica fora da media e das correlacoes, e tambem uma situacao em SITUACOES_RUPTURA (padrao:
# sem_equilibrio_geostatico e fs_abaixo_de_F0); outra situacao != ok (o calculo nao terminou) fica de fora. LIMITE_FS
# (padrao 1) troca o limiar de P_f.
export LC_ALL=C # numeros com ponto decimal em qualquer locale (awk, printf)
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUP=",${SITUACOES_RUPTURA:-sem_equilibrio_geostatico,fs_abaixo_de_F0},"
[ $# -gt 0 ] || { sed -n '2,/^export/{/^#/p}' "$0"; exit 2; }
for D in "$@"; do
    "$PKG/junta_mc.sh" "$D" > /dev/null
    [ -s "$D/mc_resumo.csv" ] || { echo "$D: nenhuma realizacao concluida" >&2; continue; }
    echo "== $(basename "$D") ($(($(wc -l < "$D/mc_resumo.csv") - 1)) realizacoes)"
    # por instante, so com mais de um instante por realizacao
    if [ -s "$D/mc_fs.csv" ] &&
        [ "$(awk -F, 'NR > 1 { n[$1]++; if (n[$1] > m) m = n[$1] } END { print m + 0 }' "$D/mc_fs.csv")" -gt 1 ]; then
        awk -F, -v L="${LIMITE_FS:-1}" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; if (!("t_h" in c) || !("FS_SRM" in c)) exit; next }
            { t = $(c["t_h"]) + 0; F = $(c["FS_SRM"]) + 0; n[t]++; if (F < L) f[t]++; if (F > 0) { m[t]++; s[t] += F; q[t] += F * F }; T[t] = 1 }
            END {
                k = 0; for (t in T) ts[++k] = t + 0
                if (!k) exit
                printf "%7s %5s %9s %8s %7s %8s %9s\n", "t (h)", "N", "media FS", "desvio", "CoV", "P_f", "CoV(P_f)"
                for (i = 1; i <= k; i++) for (j = i + 1; j <= k; j++) if (ts[j] < ts[i]) { a = ts[i]; ts[i] = ts[j]; ts[j] = a }
                for (i = 1; i <= k; i++) {
                    t = ts[i]; mu = m[t] ? s[t] / m[t] : 0; dp = m[t] > 1 ? sqrt((q[t] - m[t] * mu * mu) / (m[t] - 1)) : 0
                    pf = f[t] / n[t]; cv = pf > 0 ? sqrt((1 - pf) / (n[t] * pf)) : -1
                    printf "%7g %5d %9.4f %8.4f %7.3f %8.4f %9s\n", t, n[t], mu, dp, mu ? dp / mu : 0, pf, cv < 0 ? "-" : sprintf("%.3f", cv)
                }
            }' "$D/mc_fs.csv"
    fi
    awk -F, -v L="${LIMITE_FS:-1}" -v rup="$RUP" 'NR == 1 {
            for (i = 1; i <= NF; i++) c[$i] = i
            nv = split("c_med phi_med k_mult_med banda_c banda_phi_deg banda_k_mult banda_yc banda_area_m2 massa_c " \
                       "massa_phi_deg massa_k_mult massa_area_m2", vars, " ")
            if (!("FS_min" in c)) { print "sem a coluna FS_min em mc_resumo.csv"; sem = 1; exit }
            next
        }
        {
            if ($(c["FS_min"]) == "") next
            sit = ("situacao" in c) ? $(c["situacao"]) : "ok"
            r = (sit != "" && index(rup, "," sit ",") > 0)
            if (sit != "" && sit != "ok" && !r) { nok++; next } # o calculo nao terminou: o FS_min nao e o da realizacao
            F = $(c["FS_min"]) + 0; if (r && F > 0) F = -1
            N++; if (F < L) f++
            if ("seg" in c) seg += $(c["seg"])
            if ("cortes" in c) cortes += $(c["cortes"])
            if ("retomadas" in c) ret += $(c["retomadas"])
            if (F <= 0) { semeq++; next }
            m++; s += F; q += F * F
            for (j = 1; j <= nv; j++) {
                v = vars[j]
                if (!(v in c) || $(c[v]) == "") continue
                x = $(c[v]) + 0; k[v]++; sx[v] += x; sy[v] += F; sxx[v] += x * x; syy[v] += F * F; sxy[v] += x * F
            }
        }
        END {
            if (sem) exit 1
            mu = m ? s / m : 0; dp = m > 1 ? sqrt((q - m * mu * mu) / (m - 1)) : 0; pf = N ? f / N : 0
            printf "FS minimo de cada realizacao: media %.4f, desvio %.4f; P_f = P(FS_min < %g) = %.4f", mu, dp, L, pf
            if (pf > 0) printf " (CoV %.3f)", sqrt((1 - pf) / (N * pf))
            printf "\n%d realizacoes (%d sem equilibrio, contadas como ruptura; mais %d com o calculo incompleto, fora desta conta)", N, semeq, nok
            if ("cortes" in c) printf ", %d cortes de passo", cortes
            if ("retomadas" in c) printf ", %d retomadas", ret
            if ("seg" in c) printf "; %.1f min por realizacao", N ? seg / N / 60 : 0
            printf "\n"
            linha = ""
            for (j = 1; j <= nv; j++) {
                v = vars[j]; if (k[v] < 3) continue
                den = (k[v] * sxx[v] - sx[v] * sx[v]) * (k[v] * syy[v] - sy[v] * sy[v])
                if (den > 0) linha = linha sprintf(" %s %.2f;", v, (k[v] * sxy[v] - sx[v] * sy[v]) / sqrt(den))
            }
            if (linha != "") printf "correlacao do FS minimo com:%s\n", linha
        }' "$D/mc_resumo.csv"
done
