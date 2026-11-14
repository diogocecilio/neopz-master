#!/bin/bash
# Teste do ../confiabilidade.wl no Mathics, o interpretador livre da linguagem do Mathematica (pip install mathics3),
# para quem nao tem o Mathematica: substitutos das funcoes estatisticas (shims.m), as definicoes do .wl (secoes 1 a 3,
# 4 ate o FORM, 5 sem os cenarios e 6) e os testes teste_p1.m e teste_p2.m. Cerca de 1 min.
set -e
cd "$(dirname "$0")"
python3 - <<'PY'
s = open('../confiabilidade.wl').read().replace('ClearAll["Global`*"];', '')
i4, i5 = s.index('(*4. Problema 1'), s.index('(*5. Problema 2')
imc, icen = s.index('mc1 = monteCarlo'), s.index('(* um cenario: referencia')
i6, i7 = s.index('(*6. As funcoes do codigo original'), s.index('(*7. Esboco')
open('defs.wl', 'w').write(s[:imc] + s[i5:icen] + s[i6:i7].rsplit('(* ::Section:: *)', 1)[0])
PY
mathics -q -c '$IterationLimit = 10^7; Get["shims.m"]; Get["defs.wl"]; Get["teste_p1.m"]; Get["teste_p2.m"]' < /dev/null
