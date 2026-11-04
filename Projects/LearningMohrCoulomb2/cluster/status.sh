#!/bin/bash
# status.sh: situacao e tempos de cada caso de casos.txt, lidos dos arquivos que o calculo grava enquanto roda
# (balanco.csv a cada passo do transiente, fs.csv a cada FS). Funciona no cluster e na sua maquina; com squeue, mostra a fila.
#   ./status.sh [pasta]           tabela (pasta de resultados; padrao ./resultados)
#   ./status.sh --csv [pasta]     a mesma informacao em CSV, para guardar
#   ./status.sh --slurm [pasta]   acrescenta o sacct do ultimo job de casos submetido: estado, duracao, memoria maxima e no
#   ./status.sh --comparar A B    compara duas pastas de resultados (por exemplo a sua e a trazida do cluster): s/passo nos
#                                 passos que as duas tem, s/FS nos instantes que as duas tem e a maior diferenca de FS_SRM
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASOS="$PKG/casos.txt"
MODO=tabela
case "${1:-}" in
  --csv) MODO=csv; shift ;;
  --slurm) MODO=slurm; shift ;;
  --comparar) MODO=comparar; shift ;;
  -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
esac

# resumo <pasta do caso>: estado,passos,passos_total,t_h,seg_transiente,newton_soma,fs,fs_total,seg_fs,host,cpu
resumo() {
  local d="$1" N="" M="" ex="-" host="-" cpu="-" arqs=()
  [ -d "$d" ] || { echo "nao iniciado,0,,,0,0,0,,0,-,-"; return; }
  if [ -f "$d/plano.txt" ]; then
    N="$(sed -n 's/.*passos=\([0-9]*\).*/\1/p' "$d/plano.txt")"
    M="$(sed -n 's/.* fs=\([0-9]*\).*/\1/p' "$d/plano.txt")"
  elif [ -f "$d/caso.txt" ]; then
    N="$(tr ' ' '\n' < "$d/caso.txt" | awk -F= '$1 == "Tf_h" { tf = $2 } $1 == "Dt_h" { dt = $2 } END { if (dt > 0) printf "%d", tf / dt + 0.5 }')"
  fi
  [ -f "$d/exit.txt" ] && ex="$(sed -n 's/^exit=//p' "$d/exit.txt")"
  if [ -s "$d/tempos.csv" ] && [ "$(wc -l < "$d/tempos.csv")" -gt 1 ]; then
    host="$(tail -n 1 "$d/tempos.csv" | cut -d, -f6)"
    cpu="$(tail -n 1 "$d/tempos.csv" | cut -d, -f8-)"
  elif [ -s "$d/maquina.txt" ]; then
    host="$(sed -n 's/^host=\([^ ]*\).*/\1/p' "$d/maquina.txt")"
    cpu="$(sed -n 's/.* cpu=//p' "$d/maquina.txt" | tr ',' ' ')"
  fi
  [ -s "$d/balanco.csv" ] && arqs+=("$d/balanco.csv")
  [ -s "$d/fs.csv" ] && arqs+=("$d/fs.csv")
  awk -F, -v N="$N" -v M="$M" -v ex="$ex" -v host="$host" -v cpu="$cpu" '
    FNR == 1 { tipo = (FILENAME ~ /balanco\.csv$/) ? "b" : "f"; nc[tipo] = NF; for (i = 1; i <= NF; i++) c[tipo "," $i] = i; next }
    NF < nc[tipo] { next }   # linha incompleta: processo interrompido no meio da escrita
    tipo == "b" { np++; th = $(c["b,t_h"]); sp += $(c["b,seg_passo"]); nw += $(c["b,newton"]) }
    tipo == "f" { nf++; sf += $(c["f,seg"]) }
    END {
      if (ex == "0") est = (M != "" && nf + 0 >= M + 0) ? "concluido" : "FS"
      else if (ex == "-") est = (np > 0) ? "transiente" : "iniciando"
      else est = "falhou(exit=" ex ")"
      printf "%s,%d,%s,%s,%.3f,%d,%d,%s,%.3f,%s,%s\n", est, np, N, th, sp, nw, nf, M, sf, host, cpu
    }' ${arqs[@]+"${arqs[@]}"} < /dev/null
}

# horas/minutos: m -> "1h05"
FMT_H='function hm(m) { if (m == "") return "-"; m = int(m + 0.5); return sprintf("%dh%02d", int(m / 60), m % 60) }'

if [ "$MODO" = comparar ]; then
  A="${1:?uso: ./status.sh --comparar <pasta A> <pasta B>}"; B="${2:?uso: ./status.sh --comparar <pasta A> <pasta B>}"
  # nos e CPUs em que os casos da pasta rodaram (maquina.txt, ou a ultima fase de tempos.csv)
  maq() {
    local d r
    for d in "$1"/*/*/; do
      if [ -s "$d/maquina.txt" ]; then r="$(sed -n 's/^host=\([^ ]*\).* cpu=\(.*\)/\1 (\2)/p' "$d/maquina.txt")"
      elif [ -f "$d/tempos.csv" ] && [ "$(wc -l < "$d/tempos.csv")" -gt 1 ]; then r="$(tail -n 1 "$d/tempos.csv" | awk -F, '{ print $6 " (" $8 ")" }')"
      else continue; fi
      echo "$r"
    done | sort | uniq -c | awk '{ c = $1; $1 = ""; printf "%s%s x%d", (NR > 1 ? ";" : ""), $0, c } END { if (!NR) printf " ?" }'
  }
  # solver linear registrado em caso.txt; resultados anteriores ao registro rodaram todos com o LU skyline
  solv() { local v; v="$(cat "$1"/*/*/caso.txt 2>/dev/null | tr ' ' '\n' | sed -n 's/^solver=//p' | sort -u | paste -sd, -)"; echo "${v:-skyline (nao registrado)}"; }
  echo "A = $A:$(maq "$A"); solver $(solv "$A")"
  echo "B = $B:$(maq "$B"); solver $(solv "$B")"
  printf '%-19s %7s %9s %9s %6s %5s %8s %8s %6s %11s\n' caso passos "s/passo A" "s/passo B" "B/A" FS "s/FS A" "s/FS B" "B/A" "max|dFS|"
  while read -r caso _; do
    [ -n "$caso" ] && [ -d "$A/$caso" ] && [ -d "$B/$caso" ] || continue
    awk -v A="$A/$caso" -v B="$B/$caso" -v caso="$caso" '
      function le_bal(arq, p,   l, h, nh, i, v, cs, n) {
        if ((getline l < arq) <= 0) return 0
        nh = split(l, h, ","); for (i = 1; i <= nh; i++) if (h[i] == "seg_passo") cs = i
        while ((getline l < arq) > 0) { if (split(l, v, ",") < nh) continue; n++; S[p "," n] = v[cs] }
        close(arq); return n + 0
      }
      function le_fs(arq, p,   l, h, nh, i, v, cf, cs) {
        if ((getline l < arq) <= 0) return
        nh = split(l, h, ","); for (i = 1; i <= nh; i++) { if (h[i] == "FS_SRM") cf = i; if (h[i] == "seg") cs = i }
        while ((getline l < arq) > 0) { if (split(l, v, ",") < nh) continue; F[p "," v[1]] = v[cf]; G[p "," v[1]] = v[cs]; T[v[1]] = 1 }
        close(arq)
      }
      function r(a, b) { return a > 0 ? sprintf("%.2f", b / a) : "-" }
      BEGIN {
        na = le_bal(A "/balanco.csv", "a"); nb = le_bal(B "/balanco.csv", "b"); n = na < nb ? na : nb
        for (i = 1; i <= n; i++) { ta += S["a," i]; tb += S["b," i] }
        le_fs(A "/fs.csv", "a"); le_fs(B "/fs.csv", "b")
        for (t in T) if ((("a," t) in F) && (("b," t) in F)) {
          m++; ga += G["a," t]; gb += G["b," t]
          d = F["a," t] - F["b," t]; if (d < 0) d = -d; if (d > dm) dm = d
        }
        printf "%-19s %7d %9s %9s %6s %5d %8s %8s %6s %11s\n", caso, n, n ? sprintf("%.3f", ta / n) : "-", n ? sprintf("%.3f", tb / n) : "-",
               r(ta, tb), m, m ? sprintf("%.1f", ga / m) : "-", m ? sprintf("%.1f", gb / m) : "-", r(ga, gb), m ? sprintf("%.2e", dm) : "-"
      }'
  done < "$CASOS"
  exit 0
fi

SAIDA="${1:-$PKG/resultados}"
# tarefas do array na fila: indice -> estado (squeue -r lista uma linha por tarefa)
declare -A FILA=() VIVO=()
EU="${USER:-$(id -un)}"
SQ=0; VIVOS=" "
if command -v squeue >/dev/null 2>&1 && squeue -h -u "$EU" -o %i >/dev/null 2>&1; then
  SQ=1
  # array (job_slurm.sbatch): uma linha por tarefa (%K = indice do caso)
  while read -r k e; do case "$k" in ''|*[!0-9]*) ;; *) FILA[$k]="$e" ;; esac; done \
    < <(squeue -h -r -u "$EU" -n rahardjo45 -o "%K %T" 2>/dev/null)
  VIVOS=" $(squeue -h -u "$EU" -n rahardjo45 -t RUNNING -o %i 2>/dev/null | tr '\n' ' ')"
  VIVOS="$VIVOS$(squeue -h -r -u "$EU" -n rahardjo45 -t RUNNING -o %i 2>/dev/null | tr '\n' ' ') "
fi
# caso em execucao: <pasta do caso>/ativo ("job no pid inicio"; ou logs/casos/<n>.ativo, versao anterior) de um job vivo, ou de
# um processo vivo neste no quando roda sem SLURM. VIVO: sim, nao (sabe-se que parou) ou ? (versao antiga, sem como saber)
EU_NO="$(hostname -s)"
k=0
while read -r caso _; do
  k=$((k + 1)); [ -n "$caso" ] || continue
  v="?"; [ "$SQ" = 1 ] && v=nao
  [ "${FILA[$k]:-}" = RUNNING ] && v=sim
  [ "$v" != sim ] && [ -f "$SAIDA/$caso/maquina.txt" ] && [ ! -f "$SAIDA/$caso/ativo" ] && v=nao
  for f in "$SAIDA/$caso/ativo" "$PKG/logs/casos/$k.ativo"; do
    [ "$v" != sim ] && [ -f "$f" ] || continue
    read -r j no pid _ < "$f"
    if [ "$SQ" = 1 ] && [[ "$VIVOS" == *" $j "* ]]; then v=sim
    elif [ "$j" = local ] && [ "$no" = "$EU_NO" ] && [[ "$pid" =~ ^[0-9]+$ ]]; then kill -0 "$pid" 2>/dev/null && v=sim || v=nao
    fi
  done
  VIVO[$k]="$v"; [ "$v" = sim ] && FILA[$k]=RUNNING
done < "$CASOS"
if [ "$SQ" = 1 ]; then
  if [ "$MODO" != csv ]; then
    squeue -h -u "$EU" -n rahardjo45_build,rahardjo45 -o "%j: job %i %T ha %M, %C CPUs, %R" 2>/dev/null \
      | sed 's/^rahardjo45_build:/compilacao:/; s/^rahardjo45:/casos:/'
    [ -x "$PKG/bin/rahardjo" ] && echo "bin/rahardjo: $(date -r "$PKG/bin/rahardjo" '+%F %T')" || echo "bin/rahardjo: ainda nao compilado"
  fi
fi

# compilacao em segundo plano no no de login (submeter.sh com COMPILAR=login)
if [ "$MODO" != csv ] && [ -f "$PKG/logs/compilar_login.pid" ]; then
  read -r pid no ini < "$PKG/logs/compilar_login.pid"
  if [ "$no" = "$(hostname -s)" ] && kill -0 "$pid" 2>/dev/null; then
    echo "compilacao no no de login ($no): rodando desde $ini $(grep -o '\[ *[0-9]*%\]' "$PKG/logs/compilar_login.log" 2>/dev/null | tail -n 1)"
  elif [ ! -x "$PKG/bin/rahardjo" ]; then
    echo "compilacao no no de login: terminou sem gerar bin/rahardjo; veja tail -30 logs/compilar_login.log"
  fi
fi

n=0
{
  while read -r caso _; do
    n=$((n + 1)); [ -n "$caso" ] || continue
    echo "$n,$caso,${FILA[$n]:--},${VIVO[$n]:-?},$(resumo "$SAIDA/$caso")"
  done < "$CASOS"
} | awk -F, -v modo="$MODO" -v saida="$SAIDA" "$FMT_H"'
  {
    n = $1; caso = $2; fila = $3; vivo = $4; est = $5; np = $6; N = $7; th = $8; sp = $9; nw = $10; nf = $11; M = $12; sf = $13
    host = $14; cpu = $15; for (i = 16; i <= NF; i++) cpu = cpu "," $i
    # sem processo vivo, um caso no transiente ou nos FS esta parado (job cancelado, limite de tempo, maquina desligada)
    if (vivo == "nao" && (est == "transiente" || est == "FS" || est == "iniciando")) est = "parado"
    spp = np ? sp / np : ""; sfs = nf ? sf / nf : ""; nwm = np ? nw / np : ""
    dec = (np || nf) ? (sp + sf) / 60 : ""
    # o que falta: passos restantes no ritmo atual + FS restantes no ritmo atual ("+?" se ainda nao houve FS)
    # sem plano.txt (resultados de versoes antigas) o total de FS e desconhecido: "?" em vez de zero
    falta = ""; mais = ""
    if (est == "concluido") falta = 0
    else if (est !~ /^falhou/ && est != "nao iniciado" && est != "parado") {
      f = 0; ok = 1
      if (est != "FS") { if (N != "" && np) f += (N - np) * spp; else ok = 0 }
      if (M == "") mais = "+?"
      else if (M - nf > 0) { if (nf) f += (M - nf) * sfs; else mais = "+?" }
      if (ok) falta = f / 60
    }
    faltatxt = falta == "" ? "-" : (falta == 0 && mais != "" ? "?" : hm(falta) mais)
    if (modo == "csv") {
      if (NR == 1) print "n,caso,fila,estado,passos,passos_total,t_h,seg_transiente,seg_por_passo,newton_medio,fs,fs_total,seg_fs,seg_por_fs,decorrido_min,falta_min,host,cpu"
      printf "%d,%s,%s,%s,%d,%s,%s,%.1f,%s,%s,%d,%s,%.1f,%s,%s,%s,%s,%s\n", n, caso, fila, est, np, N, th, sp,
             np ? sprintf("%.4f", spp) : "", np ? sprintf("%.2f", nwm) : "", nf, M, sf, nf ? sprintf("%.2f", sfs) : "",
             dec == "" ? "" : sprintf("%.1f", dec), falta == "" ? "" : sprintf("%.1f", falta) mais, host, cpu
      next
    }
    if (NR == 1) printf "%2s %-19s %-15s %-8s %9s %5s %7s %6s %7s %6s %9s %8s %s\n", "n", "caso", "estado", "fila", "passo", "t(h)",
                        "s/passo", "newton", "FS", "s/FS", "decorrido", "falta", "maquina"
    printf "%2d %-19s %-15s %-8s %9s %5s %7s %6s %7s %6s %9s %8s %s\n", n, caso, est, fila, np "/" (N == "" ? "?" : N), th == "" ? "-" : sprintf("%.1f", th),
           np ? sprintf("%.3f", spp) : "-", np ? sprintf("%.1f", nwm) : "-", nf "/" (M == "" ? "?" : M), nf ? sprintf("%.1f", sfs) : "-",
           hm(dec), faltatxt, host
    conta[est ~ /^falhou/ ? "falhou" : est]++; tot += dec; if (host != "-") maq[host] = cpu
  }
  END {
    if (modo == "csv") exit
    printf "\nconcluidos %d, em FS %d, no transiente %d, iniciando %d, parados %d, falharam %d, nao iniciados %d de %d; CPU usada %s\n",
           conta["concluido"], conta["FS"], conta["transiente"], conta["iniciando"], conta["parado"], conta["falhou"], conta["nao iniciado"],
           NR, hm(tot)
    for (h in maq) printf "  %s: %s\n", h, maq[h]
    printf "(%s; passo = passos feitos/total do transiente; FS = instantes feitos/total)\n", saida
  }'

if [ "$MODO" = slurm ] && command -v sacct >/dev/null 2>&1 && [ -s "$PKG/logs/jobs.txt" ]; then
  A="$(tail -n 1 "$PKG/logs/jobs.txt" | sed -n 's/.*casos=\([0-9]*\).*/\1/p')"
  echo; echo "sacct do job $A (job ou tarefa, estado, duracao, memoria maxima, no):"
  if ! saida="$(sacct -j "$A" -n -P -o JobID,State,Elapsed,MaxRSS,NodeList 2>&1)" || [[ "$saida" == *error* ]]; then
    echo "  sacct indisponivel: $(echo "$saida" | tail -n 1)"
  else
    echo "$saida" | awk -F'|' '$1 ~ /\.batch$/ { sub(/\.batch$/, "", $1); printf "  %-14s %-11s %-10s %-9s %s\n", $1, $2, $3, $4, $5 }'
  fi
fi
