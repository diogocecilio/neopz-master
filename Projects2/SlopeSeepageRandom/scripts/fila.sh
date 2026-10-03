#!/bin/bash
# Executa os comandos de um arquivo (um por linha; linhas vazias e # são ignoradas) com P processos em paralelo.
# Uso: fila.sh <arquivo de comandos> <processos>
set -u
grep -v '^\s*\(#\|$\)' "$1" | xargs -P "$2" -I{} bash -c '{}'
