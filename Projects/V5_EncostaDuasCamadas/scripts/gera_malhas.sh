#!/bin/bash
# Gera as malhas do V5 a partir de malhas/encosta.geo (gmsh 4.x; GMSH=<executavel> se nao estiver no PATH):
#   encosta.msh        a malha do caso completo (16 linhas na camada superior, 0,25 m na zona; ~4300 nos)
#   encosta_grossa.msh a malha dos testes rapidos (6 linhas, 0,6 m; ~1000 nos)
#   encosta_muito_fina.msh estudo de malha (24 linhas, 0,15 m)
set -euo pipefail
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")/../malhas" && pwd)"
GMSH="${GMSH:-gmsh}"
cd "$AQUI"
"$GMSH" -2 encosta.geo -format msh41 -o encosta.msh > /dev/null
"$GMSH" -2 encosta.geo -format msh41 -setnumber nz 6 -setnumber hx 0.6 -setnumber hxe 1.5 -setnumber hxd 1.0 \
    -setnumber hz 0.7 -setnumber hb 4 -o encosta_grossa.msh > /dev/null
"$GMSH" -2 encosta.geo -format msh41 -setnumber nz 24 -setnumber hx 0.15 -setnumber hxe 0.5 -setnumber hxd 0.3 \
    -setnumber hz 0.2 -setnumber hb 2.5 -o encosta_muito_fina.msh > /dev/null
ls -la "$AQUI"/*.msh
