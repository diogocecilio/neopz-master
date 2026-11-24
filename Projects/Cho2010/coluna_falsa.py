#!/usr/bin/python3 -IS
"""coluna_falsa.py: executavel falso com a interface do SlopeRainColumn que o line_sampling.py usa (prog=coluna,
falha=fs), para testar o LS em estagios em segundos (teste T0 da nota viabilidade_ls_adaptativo.md).

    coluna_falsa.py n=1 primeira=0 xi=<arq> M=<K> alvo=<arq> [beta=2.2] [escala=0.1] ... tempos=<T> saida=<rotulo>

Le os xi (um por linha, como o line_sampling.py grava) e a direcao a do alvo= (um numero por linha, |a| = 1), e grava
<rotulo>_0_0.csv (t_h, situacao, t_passagem, FS_passagem, t_parou) com uma linha em t_h = T (o ultimo tempos=) e
    FS = exp(escala (beta - a.xi))
(estado limite linear: falha em a.xi > beta; P_f exata = Phi(-beta)). As opcoes alvo=, beta= e escala= vao na lista
"exe" do form.json. Sem imports alem de math e sys (partida rapida).
"""
import math
import sys

op = {}
for a in sys.argv[1:]:
    k, v = a.split("=", 1)
    op[k] = v  # a ultima vence (o line_sampling.py poe tempos=T depois das opcoes do form.json)
x = [float(t) for t in open(op["xi"]).read().split()]
a = [float(t) for t in open(op["alvo"]).read().split()]
g = float(op.get("beta", 2.2)) - sum(p * q for p, q in zip(a, x))
fs = math.exp(float(op.get("escala", 0.1)) * g)
T = float(op["tempos"].split(",")[-1])
with open(op["saida"] + "_0_0.csv", "w") as f:
    f.write("t_h,situacao,t_passagem,FS_passagem,t_parou\n")
    f.write(f"{T!r},ok,-1,{fs!r},-1\n")
