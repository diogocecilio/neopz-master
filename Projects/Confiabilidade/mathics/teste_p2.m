(* problema 2: Nataf, FORM com e sem Armijo, G vetorizada e amostragem por importancia no cenario A *)
Print["rhoZ LN-LN: ", rhoZ[lnDist[10., 0.3], lnDist[30., 0.1], -0.5], " (Python: -0.516123)"];
Do[Module[{m = mod2[mu], G, f, fp},
   G[u_?VectorQ] := g2X[xDeU[m, u]];
   f = form[G, 3]; fp = form[G, 3, "Armijo" -> False];
   Print["mu_h = ", mu, ": Armijo beta ", f["beta"], " (", f["iteracoes"], " it, conv ", f["convergiu"], ") alfa^2 ",
     f["alfa"]^2, " | puro beta ", fp["beta"], " (", fp["iteracoes"], " it, conv ", fp["convergiu"], ")"]],
  {mu, {-3., -2., -0.5}}];
(* Python: -3: 5.19746 (12 it), puro 4.58648 sem convergir; -2: 3.54141 (21 it), puro 35 it; -0.5: 1.07813 (7 it) *)
Module[{m = mod2[-2.], G, Gv, f, r},
  G[u_?VectorQ] := g2X[xDeU[m, u]];
  Gv[us_?MatrixQ] := g2X[Transpose[xDeUs[m, us]]];
  f = form[G, 3];
  Print["Gv contra G: ", Gv[{{0.1, -0.2, 0.3}, {1., 1., 1.}}], " ", {G[{0.1, -0.2, 0.3}], G[{1., 1., 1.}]}];
  r = importancia[Gv, f["u"], 400];
  Print["importancia, cenario A, N = 400: pf = ", r["pf"], " (CoV ", r["cov"], "); referencia 1.535*10^-4"]];
