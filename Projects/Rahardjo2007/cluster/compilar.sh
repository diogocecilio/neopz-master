#!/bin/bash
# compilar.sh: compila so o Rahardjo2007 numa build e num bin proprios (build/ e bin/ desta pasta), com o compilar.sh do
# SlopeRandomFieldRefinement/cluster (modulos, MKL sequencial, versao em bin/versao.txt). No no de login do fermi (a
# primeira vez leva uns 15 a 30 min):  ./compilar.sh
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="${BUILD:-$AQUI/build}" BIN="$AQUI/bin" ALVOS=Rahardjo2007/Rahardjo2007 NOMES=rah2007 \
    exec "$AQUI/../../SlopeRandomFieldRefinement/cluster/compilar.sh"
