(* ===== KL ANALITICO (livro, eq. 5.54-5.60), dominio [-a,a]^2 com a = 0.5, c = 1/b = 1 ===== *)
a = 0.5; c = 1./LX;
nw = 10;   (* raizes 1D *)
(* eq. 5.59: i impar: c - w tan(w a) = 0 ; i par: w + c tan(w a) = 0. A i-esima raiz fica em ((i-1) pi/(2a), i pi/(2a)) *)
om = Table[If[OddQ[i],
     x /. FindRoot[c - x Tan[x a], {x, (i - 1) Pi/(2 a) + 0.01, (i - 1) Pi/(2 a) + 0.01, i Pi/(2 a) - 0.01}],
     x /. FindRoot[x + c Tan[x a], {x, (i - 0.5) Pi/(2 a) + 0.01, (i - 1) Pi/(2 a) + 0.01, i Pi/(2 a) - 0.01}]], {i, nw}];
lam1D = 2 c/(om^2 + c^2);                                              (* eq. 5.56 *)
f1D[i_, x_] := If[OddQ[i], Cos[om[[i]] x]/Sqrt[a + Sin[2 om[[i]] a]/(2 om[[i]])],   (* eq. 5.57 *)
                           Sin[om[[i]] x]/Sqrt[a - Sin[2 om[[i]] a]/(2 om[[i]])]];  (* eq. 5.58 *)
(* modos 2D: f_i(x) f_j(y), lambda = lambda_i lambda_j (eq. 5.50-5.51), ordenados *)
pairs = SortBy[Flatten[Table[{i, j}, {i, nw}, {j, nw}], 1], -lam1D[[#[[1]]]] lam1D[[#[[2]]]] &];
lamA = (lam1D[[#[[1]]]] lam1D[[#[[2]]]]) & /@ pairs;
phiA = Table[f1D[p[[1]], #[[1]]] f1D[p[[2]], #[[2]]] & /@ coords, {p, Take[pairs, sz]}];

(* conferencias *)
Print["1a raiz w1 = ", om[[1]], "   lambda 1D = ", Take[lam1D, 4]];
Print["checagem integral 1D: ", NIntegrate[Exp[-c Abs[0.2 - y]] f1D[1, y], {y, -a, a}]/(lam1D[[1]] f1D[1, 0.2])];
TableForm[Transpose[{Range[8], Take[lamA, 8], Take[lam, 8]}],
  TableHeadings -> {None, {"modo", "analitico", "Nystrom"}}] // Print;
Print["variancia de E no canto A (deveria ser 1) M=2,4,25: analitico ",
  Table[Sum[lamA[[i]] phiA[[i, nA]]^2, {i, m}], {m, {2, 4, sz}}], "  Nystrom ",
  Table[Sum[lam[[i]] phi[[i, nA]]^2, {i, m}], {m, {2, 4, sz}}]];

(* matrizes K^(i) com as autofuncoes analiticas e Monte Carlo com os MESMOS xi *)
KkA = Table[First@AssembleE[mesh, zero, Sqrt[lamA[[i]]] phiA[[i]]], {i, sz}];
KkAf = Normal[#][[free, free]] & /@ KkA;
MCrunA[sig_?NumericQ, M_Integer] := Module[{u},
   u = Table[LinearSolve[K0f + sig Total[xi[[s, 1 ;; M]] KkAf[[1 ;; M]]], Ff][[posA]], {s, nMC}];
   {Mean[u], StandardDeviation[u]}];
resA = Table[{s, MCrunA[s, 2][[2]], MCrunA[s, 4][[2]], MCrunA[s, sz][[2]]}, {s, sigs}];

(* o livro plota sigma_u NORMALIZADO por sigma_max = 0.433 (legenda da Fig. 5.46) *)
bookAbs = {#[[1]], 0.433 #[[2]]} & /@ book;
tab = Table[{sigs[[i]], res[[i, 2]], res[[i, 4]], resA[[i, 2]], resA[[i, 4]],
    If[sigs[[i]] >= 0.1, 0.433 Interpolation[book, InterpolationOrder -> 1][sigs[[i]]], "-"]}, {i, Length[sigs]}];
Print[TableForm[Map[If[NumberQ[#], NumberForm[#, {4, 3}], #] &, tab, {2}],
   TableHeadings -> {None, {"sigma_E", "Nys M=2", "Nys M=25", "Anal M=2", "Anal M=25", "livro x0.433"}}]];

gAbs = Show[ListLinePlot[{res[[All, {1, 2}]], res[[All, {1, 4}]], resA[[All, {1, 2}]], resA[[All, {1, 4}]]},
    PlotMarkers -> Automatic, PlotStyle -> {Blue, Darker@Green, {Blue, Dashed}, {Darker@Green, Dashed}},
    PlotLegends -> {"Nystrom M=2", "Nystrom M=25", "KL analitico M=2", "KL analitico M=25"}],
   ListPlot[bookAbs, PlotStyle -> Black, PlotMarkers -> {"*", 18}, PlotLegends -> {"Livro MCS \[Times] 0.433"}],
   Frame -> True, FrameLabel -> {"desvio padrao de E", "desvio padrao de u_y(A)"},
   PlotRange -> {{0, 0.31}, {0, 0.5}}, GridLines -> Automatic, ImageSize -> 600]
