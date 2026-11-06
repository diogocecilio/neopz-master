# Ajuste para o seu cluster (veja `module avail`). Lido por submeter.sh, pela compilacao e pelos jobs.
# fermi (CESUP/UFRGS): gcc do sistema + MKL do oneAPI 2025.1 (o modulo define MKLROOT, que o FindMKL usa).

# Em shells nao interativos o comando module pode nao estar definido
if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash /etc/profile.d/modules.sh; do
    [ -r "$f" ] && { . "$f"; break; }
  done
fi

# carregar <modulo>: os modulos do oneAPI declaram pre-requisitos (prereq) que o Lmod nao carrega sozinho; ele so
# avisa "At least one of these module(s) must be loaded: X". Aqui X e carregado antes e o modulo e tentado de novo.
carregar() {
  local pilha="$1" alvo dep err n=0
  err="$(mktemp)"
  while [ -n "$pilha" ] && [ "$n" -lt 12 ]; do
    n=$((n + 1))
    alvo="${pilha%% *}"
    if module load "$alvo" 2>"$err"; then
      pilha="${pilha#"$alvo"}"; pilha="${pilha# }"
    else
      dep="$(awk '/must be loaded/ { getline; print $1; exit }' "$err")"
      [ -n "$dep" ] || break
      pilha="$dep $pilha"
    fi
  done
  [ -z "$pilha" ] || cat "$err" >&2
  rm -f "$err"
  [ -z "$pilha" ]
}

carregar mkl/2025.1
# fermi: os nos de computacao nao tem cmake (o de login tem; c++ e gcc 14 no login, gcc 15 nos nos). submeter.sh compila
# no no de login, em segundo plano, e submete os casos quando terminar. Em outro cluster: COMPILAR=job (compila num job).
export COMPILAR="${COMPILAR:-login}"
# export USAR_MKL=OFF      # sem MKL: descomente e comente a linha acima (skyline, mais lento)
# export PARTICAO=normal   # particao do SLURM; vazio usa a padrao
