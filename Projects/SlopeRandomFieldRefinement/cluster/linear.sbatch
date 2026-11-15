#!/bin/bash
# linear.sbatch: teste da amostragem por importancia SEM elementos finitos (etapa 2 do plano de trabalho), num no.
# O Cho2010 com linear=beta alfa=<direcao> usa o sorteio real dos campos (semente, desloc_xi, peso), mas a falha e
# alfa . xi > beta, com P_f exata = Phi(-beta). Direcao: o ponto de projeto do FORM M = 20 da campanha etapa1. Para beta =
# 2, 3 e 4: Monte Carlo direto (10^6) e AI (10^4, semente 3) com o deslocamento do FORM M = 20 (centrado na direcao
# certa) e do M = 10 (um pouco fora). Resultado em $DATA/etapa2_linear/resultado.txt.
#   sbatch linear.sbatch
#SBATCH --job-name=etapa1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=logs/%x_linear_%j.out
[ -r /etc/profile.d/scratch.sh ] && source /etc/profile.d/scratch.sh
cd "$SLURM_SUBMIT_DIR" || exit 1
source ../../LearningMohrCoulomb2/cluster/modulos.sh
X="$PWD/bin/Cho2010" P="$PWD/../posproc/amostragem_importancia.py"
E="${DATA:-$HOME}/etapa1" S="${DATA:-$HOME}/etapa2_linear"
mkdir -p "$S" && cd "$S" || exit 1
export PZ_THREADS=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for m in 10 20; do
    [ -s "$E/form_M$m/ustar.txt" ] || { echo "falta $E/form_M$m/ustar.txt" >&2; exit 1; }
    python3 -I - "$E/form_M$m/ustar.txt" "$m" <<'PY'
import math, sys
u = [float(x) for x in open(sys.argv[1]).read().split()]
m = sys.argv[2]
n = math.sqrt(sum(x * x for x in u))
open(f"alfa_M{m}.txt", "w").write("".join(repr(x / n) + "\n" for x in u))
for b in (2, 3, 4):
    open(f"mu_M{m}_b{b}.txt", "w").write("".join(repr(b * x / n) + "\n" for x in u))
PY
done
for b in 2 3 4; do
    mkdir -p b$b/direto b$b/ai20 b$b/ai10
    "$X" n=1000000 linear=$b alfa=alfa_M20.txt saida=b$b/direto/mc > b$b/direto.log 2>&1 &
    "$X" n=10000 linear=$b alfa=alfa_M20.txt desloc_xi=mu_M20_b$b.txt semente=3 saida=b$b/ai20/mc > b$b/ai20.log 2>&1 &
    "$X" n=10000 linear=$b alfa=alfa_M20.txt desloc_xi=mu_M10_b$b.txt semente=3 saida=b$b/ai10/mc > b$b/ai10.log 2>&1 &
done
wait
{ for b in 2 3 4; do
      python3 -I -c "import math; print('== beta = $b: P_f exata = %.4e' % (0.5 * math.erfc($b / math.sqrt(2))))"
      echo "  direto (10^6): $(grep 'P(FS < 1)' b$b/direto.log | sed 's/^ *//')"
      for a in ai20 ai10; do
          echo "  $a (10^4):"
          python3 -I -B "$P" estima b$b/$a 2>&1 | grep -E 'deslocamento|AI:|P_f \(' | sed 's/^/    /'
      done
  done; } > resultado.txt
cat resultado.txt
