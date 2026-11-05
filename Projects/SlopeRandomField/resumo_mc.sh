#!/bin/bash
# resumo_mc.sh <pasta da combinacao> [...]: estatistica do Monte Carlo de cada combinacao (lida de mc_fs.csv e
# mc_resumo.csv): por instante, numero de realizacoes, media, desvio e CoV do FS, probabilidade de ruptura
# P_f = P(FS < 1) e o CoV da estimativa, sqrt((1 - P_f) / (N P_f)) (Vargas Ceron et al. 2025, eq. 68); depois o mesmo
# para o FS minimo de cada realizacao. FS = -1 (sem equilibrio nem com c'/0.6) conta como ruptura e fica fora da media.
for D in "$@"; do
    [ -s "$D/mc_fs.csv" ] || { echo "$D: sem mc_fs.csv" >&2; continue; }
    echo "== $(basename "$D") ($(($(wc -l < "$D/mc_resumo.csv") - 1)) realizacoes)"
    awk -F, 'NR > 1 { t = $2 + 0; n[t]++; if ($3 < 1) f[t]++; if ($3 > 0) { m[t]++; s[t] += $3; q[t] += $3 * $3 } ; T[t] = 1 }
        END {
            printf "%7s %5s %9s %8s %7s %8s %9s\n", "t (h)", "N", "media FS", "desvio", "CoV", "P_f", "CoV(P_f)"
            k = 0; for (t in T) ts[++k] = t + 0
            for (i = 1; i <= k; i++) for (j = i + 1; j <= k; j++) if (ts[j] < ts[i]) { a = ts[i]; ts[i] = ts[j]; ts[j] = a }
            for (i = 1; i <= k; i++) {
                t = ts[i]; mu = m[t] ? s[t] / m[t] : 0; dp = m[t] > 1 ? sqrt((q[t] - m[t] * mu * mu) / (m[t] - 1)) : 0
                pf = f[t] / n[t]; cv = pf > 0 ? sqrt((1 - pf) / (n[t] * pf)) : -1
                printf "%7g %5d %9.4f %8.4f %7.3f %8.4f %9s\n", t, n[t], mu, dp, mu ? dp / mu : 0, pf, cv < 0 ? "-" : sprintf("%.3f", cv)
            }
        }' "$D/mc_fs.csv"
    awk -F, 'NR > 1 { N++; if ($2 < 1) f++; if ($2 > 0) { m++; s += $2; q += $2 * $2 }; if ($4 != "ok") nok++; seg += $8; cortes += $6 }
        END {
            mu = m ? s / m : 0; dp = m > 1 ? sqrt((q - m * mu * mu) / (m - 1)) : 0; pf = N ? f / N : 0
            printf "FS minimo de cada realizacao: media %.4f, desvio %.4f; P_f = %.4f", mu, dp, pf
            if (pf > 0) printf " (CoV %.3f)", sqrt((1 - pf) / (N * pf))
            printf "\n%d realizacoes, %d com o transiente interrompido, %d cortes de passo; %.1f min por realizacao\n",
                   N, nok, cortes, N ? seg / N / 60 : 0
        }' "$D/mc_resumo.csv"
done
