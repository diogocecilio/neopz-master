(* ===== PROBABILIDADE DE FALHA:  g = delta - u_y(A),  falha se g <= 0 ===== *)
sigE = 0.2;          (* desvio padrao de E *)
nPf = 20000;         (* amostras de Monte Carlo *)
SeedRandom[2026];
xiPf = RandomVariate[NormalDistribution[], {nPf, sz}];

(* amostra u_y(A) e E minimo nodal para uma base KL (lambda, phi, K^(i) livres) e M termos *)
SampleU[sig_?NumericQ, M_Integer, lamB_, phiB_, KfB_] := Module[{u, Emin},
   u = Table[LinearSolve[K0f + sig Total[xiPf[[s, 1 ;; M]] KfB[[1 ;; M]]], Ff][[posA]], {s, nPf}];
   Emin = Min /@ (1. + sig xiPf[[All, 1 ;; M]] . (Sqrt[lamB[[1 ;; M]]] phiB[[1 ;; M]]));
   {u, Emin}];

(* funcao de estado limite; amostra com E <= 0 em algum no tambem conta como falha *)
gLS[u_, delta_] := delta - u;
PfMC[{u_, Emin_}, delta_] := Module[{fail, pf, n = Length[u]},
   fail = MapThread[Boole[gLS[#1, delta] <= 0 || #2 <= 0] &, {u, Emin}];
   pf = Mean[N@fail];
   {pf, If[pf > 0, Sqrt[(1 - pf)/(n pf)], Infinity]}];    (* {Pf, coef. de variacao do estimador} *)

(* amostras: Nystrom e KL analitico, M = 2 e M = 25 *)
sN2 = SampleU[sigE, 2, lam, phi, Kkf];    sN25 = SampleU[sigE, sz, lam, phi, Kkf];
sA2 = SampleU[sigE, 2, lamA, phiA, KkAf]; sA25 = SampleU[sigE, sz, lamA, phiA, KkAf];
Print["amostras com E<=0 (Nys25, Anal25): ", Count[sN25[[2]], _?(# <= 0 &)], ", ", Count[sA25[[2]], _?(# <= 0 &)]];
Print["media e desvio de u (Nys25): ", {Mean[sN25[[1]]], StandardDeviation[sN25[[1]]]}];

deltas = {1.2, 1.3, 1.4, 1.5, 1.6};
tabPf = Table[Join[{d}, Flatten[PfMC[#, d] & /@ {sN2, sN25, sA2, sA25}]], {d, deltas}];
TableForm[Map[NumberForm[#, 3] &, tabPf, {2}],
  TableHeadings -> {None, {"delta", "Pf Nys2", "cv", "Pf Nys25", "cv", "Pf Anal2", "cv", "Pf Anal25", "cv"}}]

(* curva Pf(delta) = P(u_y(A) > delta) *)
dd = Range[1.05, 1.8, 0.025];
curva[s_] := Table[{d, First@PfMC[s, d]}, {d, dd}];
gPf = ListLogPlot[{curva[sN2], curva[sN25], curva[sA2], curva[sA25]}, Joined -> True,
   PlotStyle -> {Blue, Darker@Green, {Blue, Dashed}, {Darker@Green, Dashed}},
   PlotLegends -> {"Nystrom M=2", "Nystrom M=25", "KL analitico M=2", "KL analitico M=25"},
   Frame -> True, FrameLabel -> {"\[Delta] (deslocamento limite)", "P_f = P[u_y(A) > \[Delta]]"},
   PlotLabel -> Row[{"\[Sigma]_E = ", sigE, ",  N = ", nPf}], GridLines -> Automatic, ImageSize -> 600]

(* histograma de u com o limite delta = 1.4 *)
Show[Histogram[{sN25[[1]], sA25[[1]]}, {0.02}, "PDF", ChartLegends -> {"Nystrom M=25", "KL analitico M=25"},
   ChartStyle -> {Opacity[0.5, Darker@Green], Opacity[0.5, Orange]}],
  Graphics[{Red, Thick, Line[{{1.4, 0}, {1.4, 4}}]}], Frame -> True, FrameLabel -> {"u_y(A)", "densidade"},
  PlotRange -> {{0.4, 2}, All}]
