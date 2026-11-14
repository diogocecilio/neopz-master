(* problema 1: FORM (ja impresso ao carregar), amostragem por importancia e Monte Carlo pequenos *)
r = importancia[g1Us, f1["u"], 1000];
Print["importancia, N = 1000: pf = ", r["pf"], " (CoV ", r["cov"], "); referencia 1.17688*10^-3"];
r = monteCarlo[g1Us, 3, 2000];
Print["Monte Carlo, N = 2000: pf = ", r["pf"], " (CoV ", r["cov"], ")"];
