#!/bin/bash
# resumo_mc.sh <pasta da combinacao> [...]: junta as realizacoes concluidas (junta_mc.sh) e mostra a estatistica do Monte
# Carlo: por instante, numero de realizacoes, media, desvio e CoV do FS, probabilidade de ruptura P_f = P(FS < 1) e o
# CoV da estimativa, sqrt((1 - P_f) / (N P_f)) (Vargas Ceron et al. 2025, eq. 68); depois o mesmo para o FS minimo de
# cada realizacao e a correlacao (Pearson) do FS minimo com as medias dos campos no solo e na banda de ruptura. FS = -1
# (sem equilibrio nem com c'/0.6) conta como ruptura e fica fora da media e das correlacoes.
export LC_ALL=C # numeros com ponto decimal em qualquer locale (awk, printf)
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for D in "$@"; do
    "$PKG/junta_mc.sh" "$D" > /dev/null
    [ -s "$D/mc_fs.csv" ] || { echo "$D: nenhuma realizacao concluida" >&2; continue; }
    echo "== $(basename "$D") ($(($(wc -l < "$D/mc_resumo.csv") - 1)) realizacoes)"
    awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        { t = $(c["t_h"]) + 0; F = $(c["FS_SRM"]) + 0; n[t]++; if (F < 1) f[t]++; if (F > 0) { m[t]++; s[t] += F; q[t] += F * F }; T[t] = 1 }
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
    awk -F, 'NR == 1 {
            for (i = 1; i <= NF; i++) c[$i] = i
            nv = split("c_med phi_med k_mult_med banda_c banda_phi_deg banda_k_mult banda_yc banda_area_m2", vars, " ")
            next
        }
        {
            if ($(c["FS_min"]) == "") next
            if ($(c["situacao"]) != "ok") { nok++; next } # transiente interrompido: o FS_min e o de antes da parada
            F = $(c["FS_min"]) + 0; N++; if (F < 1) f++
            seg += $(c["seg"]); cortes += $(c["cortes"]); if ("retomadas" in c) ret += $(c["retomadas"])
            if (F <= 0) next
            m++; s += F; q += F * F
            for (j = 1; j <= nv; j++) {
                v = vars[j]
                if (!(v in c) || $(c[v]) == "") continue
                x = $(c[v]) + 0; k[v]++; sx[v] += x; sy[v] += F; sxx[v] += x * x; syy[v] += F * F; sxy[v] += x * F
            }
        }
        END {
            mu = m ? s / m : 0; dp = m > 1 ? sqrt((q - m * mu * mu) / (m - 1)) : 0; pf = N ? f / N : 0
            printf "FS minimo de cada realizacao: media %.4f, desvio %.4f; P_f = %.4f", mu, dp, pf
            if (pf > 0) printf " (CoV %.3f)", sqrt((1 - pf) / (N * pf))
            printf "\n%d realizacoes (mais %d com o transiente interrompido, fora desta conta), %d cortes de passo, %d retomadas; %.1f min por realizacao\n",
                   N, nok, cortes, ret, N ? seg / N / 60 : 0
            linha = ""
            for (j = 1; j <= nv; j++) {
                v = vars[j]; if (k[v] < 3) continue
                den = (k[v] * sxx[v] - sx[v] * sx[v]) * (k[v] * syy[v] - sy[v] * sy[v])
                if (den > 0) linha = linha sprintf(" %s %.2f;", v, (k[v] * sxy[v] - sx[v] * sy[v]) / sqrt(den))
            }
            if (linha != "") printf "correlacao do FS minimo com:%s\n", linha
        }' "$D/mc_resumo.csv"
done
