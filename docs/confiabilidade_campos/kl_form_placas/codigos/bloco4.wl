(* ===== FORM (HLRF) para g(xi) = delta - u_y(A; xi) ===== *)
(* variaveis aleatorias: xi_1..xi_M ~ N(0,1) independentes (as mesmas do KL) *)

(* resolve K(xi) u = F e devolve {u_A, gradiente de u_A em xi} pelo metodo adjunto:
   K = K0 + sig Sum xi_k K_k  ->  du/dxi_k = -K^-1 (sig K_k u)
   d u_A/dxi_k = -sig lamb.(K_k u),  com K lamb = e_A  (K simetrica)  *)
SolveGrad[xiv_, sig_, KfB_] := Module[{M = Length[xiv], K, LS, u, eA, lamb},
   K = K0f + sig Total[xiv KfB[[1 ;; M]]];
   LS = LinearSolve[K];
   u = LS[Ff];
   eA = UnitVector[Length[Ff], posA];
   lamb = LS[eA];
   {u[[posA]], Table[-sig lamb . (KfB[[k]] . u), {k, M}]}];

(* roda o FORM para um delta, um sigma_E, M termos e uma base KL *)
FORMrun[delta_, sig_, M_, KfB_, iter_ : 10] := Module[{gf, gradf, res},
   gf[x_] := delta - First@SolveGrad[x, sig, KfB];        (* g = delta - u_A *)
   gradf[x_] := -Last@SolveGrad[x, sig, KfB];             (* grad g = - grad u_A *)
   res = Block[{Print}, HLRF[ConstantArray[0., M], IdentityMatrix[M], ConstantArray[0., M],
       ConstantArray["NORMAL", M], IdentityMatrix[M], gf, gradf, iter]];   (* Block[{Print}] silencia *)
   Append[res, "Pf" -> CDF[NormalDistribution[], -res["beta"]]]];

(* checagem do gradiente adjunto por diferenca finita *)
With[{x0 = RandomVariate[NormalDistribution[], 4], h = 10^-6},
  Print["grad adjunto: ", Last@SolveGrad[x0, sigE, Kkf]];
  Print["grad dif.fin: ", Table[(First@SolveGrad[x0 + h UnitVector[4, k], sigE, Kkf] -
       First@SolveGrad[x0, sigE, Kkf])/h, {k, 4}]]];

(* exemplo com as iteracoes impressas: delta = 1.4, Nystrom, M = 25 *)
r = HLRF[ConstantArray[0., sz], IdentityMatrix[sz], ConstantArray[0., sz], ConstantArray["NORMAL", sz],
   IdentityMatrix[sz], (1.4 - First@SolveGrad[#, sigE, Kkf]) &, (-Last@SolveGrad[#, sigE, Kkf]) &, 8];
Print["beta = ", r["beta"], "   Pf = ", CDF[NormalDistribution[], -r["beta"]]];
Print["alpha^2 (8 maiores) = ", Take[r["alpha2Y"], 8]];

(* tabela FORM x Monte Carlo *)
tabFORM = Table[Module[{fN2, fN25, fA2, fA25},
    fN2 = FORMrun[d, sigE, 2, Kkf]; fN25 = FORMrun[d, sigE, sz, Kkf];
    fA2 = FORMrun[d, sigE, 2, KkAf]; fA25 = FORMrun[d, sigE, sz, KkAf];
    {d, fN25["beta"], fN2["Pf"], First@PfMC[sN2, d], fN25["Pf"], First@PfMC[sN25, d],
     fA2["Pf"], First@PfMC[sA2, d], fA25["Pf"], First@PfMC[sA25, d]}], {d, deltas}];
TableForm[Map[NumberForm[#, 3] &, tabFORM, {2}],
  TableHeadings -> {None, {"delta", "beta Nys25", "FORM Nys2", "MC Nys2", "FORM Nys25", "MC Nys25",
     "FORM Anal2", "MC Anal2", "FORM Anal25", "MC Anal25"}}]

(* curva Pf(delta): FORM (linhas) x MC (pontos), Nystrom M=25 *)
ddF = Range[1.1, 1.8, 0.05];
formCurve = Table[{d, FORMrun[d, sigE, sz, Kkf]["Pf"]}, {d, ddF}];
gFORM = Show[ListLogPlot[{formCurve}, Joined -> True, PlotStyle -> {Red}, PlotLegends -> {"FORM, Nystrom M=25"}],
  ListLogPlot[{curva[sN25]}, PlotStyle -> {Darker@Green}, PlotMarkers -> {Automatic, 8},
   PlotLegends -> {"MC 20000, Nystrom M=25"}],
  Frame -> True, FrameLabel -> {"\[Delta]", "P_f"}, GridLines -> Automatic, PlotLabel -> Row[{"\[Sigma]_E = ", sigE}]]
