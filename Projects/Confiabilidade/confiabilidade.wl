(* ::Package:: *)

(* ::Title:: *)
(*Confiabilidade: FORM (HL-RF), Monte Carlo e amostragem por importancia*)


(* ::Text:: *)
(*Versao revisada do codigo original: a lista dos erros corrigidos e a teoria estao em doc/confiabilidade.pdf.*)
(*Convencoes:*)
(*  X  variaveis fisicas, com marginais F_i e correlacao RX;*)
(*  Z  normais padrao correlacionadas, z_i = Phi^-1(F_i(x_i)), com correlacao RZ (a de Nataf);*)
(*  U  normais padrao independentes, z = L.u com RZ = L.Transpose[L] (L triangular inferior); u = L^-1.z.*)
(*Falha: g(x) <= 0; no espaco U, G(u) = g(x(u)).*)
(*Os numeros de referencia nos comentarios vem da versao em Python (python/confiab.py), o mesmo algoritmo; as*)
(*estimativas de Monte Carlo e de importancia mudam com a semente, dentro do coeficiente de variacao impresso.*)


ClearAll["Global`*"];
SeedRandom[20261114];
normal01 = NormalDistribution[0, 1];


(* ::Section:: *)
(*1. Marginais e transformacao de Nataf*)


(* lognormal pela media e pelo coeficiente de variacao (o Mathematica usa os parametros de ln X) *)
lnDist[media_, cov_] := LogNormalDistribution[Log[media] - Log[1 + cov^2]/2, Sqrt[Log[1 + cov^2]]];

(* x = F^-1(Phi(z)), a transformacao exata: formas fechadas (vetorizadas) para a normal e a lognormal; para as outras
   marginais, ponto a ponto, com a cauda superior pela funcao de sobrevivencia (Phi(z) perto de 1 perde precisao) *)
xDeZ[NormalDistribution[m_, s_], z_] := m + s z;
xDeZ[LogNormalDistribution[m_, s_], z_] := Exp[m + s z];
xDeZ[d_, z_List] := xDeZ[d, #] & /@ z;
xDeZ[d_, z_?NumericQ] := If[z <= 0, InverseCDF[d, CDF[normal01, z]], InverseSurvivalFunction[d, CDF[normal01, -z]]];

zDeX[NormalDistribution[m_, s_], x_] := (x - m)/s;
zDeX[LogNormalDistribution[m_, s_], x_] := (Log[x] - m)/s;
zDeX[d_, x_?NumericQ] := InverseCDF[normal01, CDF[d, x]];

(* correlacao no espaco Z que reproduz a correlacao r das variaveis fisicas (Liu e Der Kiureghian 1986): formas fechadas
   para normal-normal, lognormal-lognormal e normal-lognormal; para outras marginais, a equacao integral *)
cvDe[d_] := StandardDeviation[d]/Mean[d];
rhoZ[_, _, 0 | 0.] := 0.;
rhoZ[_NormalDistribution, _NormalDistribution, r_] := r;
rhoZ[di_LogNormalDistribution, dj_LogNormalDistribution, r_] := Log[1 + r cvDe[di] cvDe[dj]]/(di[[2]] dj[[2]]);
rhoZ[_NormalDistribution, dj_LogNormalDistribution, r_] := r cvDe[dj]/dj[[2]];
rhoZ[di_LogNormalDistribution, _NormalDistribution, r_] := r cvDe[di]/di[[2]];
rhoZ[di_, dj_, r_] := Module[{rx, s},
   rx[s_?NumericQ] := NIntegrate[
      (xDeZ[di, z1] - Mean[di]) (xDeZ[dj, s z1 + Sqrt[1 - s^2] z2] - Mean[dj]) PDF[normal01, z1] PDF[normal01, z2],
      {z1, -8, 8}, {z2, -8, 8}, Method -> {Automatic, "SymbolicProcessing" -> 0}]/
     (StandardDeviation[di] StandardDeviation[dj]);
   s /. FindRoot[rx[s] == r, {s, r, 0.95 r}]];

(* modelo probabilistico: marginais, correlacoes e o fator de Cholesky *)
modelo[dists_List, rx_?MatrixQ] := Module[{n = Length[dists], rz, L},
   rz = Table[If[i == j, 1., rhoZ[dists[[i]], dists[[j]], rx[[i, j]]]], {i, n}, {j, n}];
   L = Transpose[CholeskyDecomposition[rz]];         (* CholeskyDecomposition da U com Transpose[U].U = rz *)
   If[Max[Abs[L . Transpose[L] - rz]] > 10^-10, Print["modelo: a decomposicao de Cholesky falhou"]];
   <|"dists" -> dists, "RX" -> rx, "RZ" -> rz, "L" -> L, "n" -> n|>];
modelo[dists_List] := modelo[dists, IdentityMatrix[Length[dists]]];

(* u -> x de um ponto; us -> xs de muitas amostras (uma por linha, vetorizado); x -> u *)
xDeU[mod_, u_?VectorQ] := MapThread[xDeZ, {mod["dists"], mod["L"] . u}];
xDeUs[mod_, us_?MatrixQ] := Transpose[MapThread[xDeZ, {mod["dists"], Transpose[us . Transpose[mod["L"]]]}]];
uDeX[mod_, x_?VectorQ] := LinearSolve[mod["L"], MapThread[zDeX, {mod["dists"], x}]];


(* ::Section:: *)
(*2. FORM: HL-RF com busca de Armijo (Zhang e Der Kiureghian 1995)*)


(* gradiente por diferencas centrais: 2n avaliacoes de G (num modelo de elementos finitos, cada uma e uma analise) *)
gradiente[G_, u_, h_] := Table[(G[u + h UnitVector[Length[u], i]] - G[u - h UnitVector[Length[u], i]])/(2 h),
   {i, Length[u]}];

(* Passo: o h das diferencas; TolG: |G|/|G(0)|; TolU: |u - (alfa.u) alfa|, u alinhado com o gradiente;
   Armijo -> False: o HL-RF puro (pode oscilar sem convergir; ver o problema 2 com mu_h <= -2,5 m) *)
Options[form] = {"Passo" -> 10.^-5, "TolG" -> 10.^-6, "TolU" -> 10.^-4, "MaxIter" -> 100, "Armijo" -> True};
form[G_, n_Integer, OptionsPattern[]] := Module[
   {h = OptionValue["Passo"], u = ConstantArray[0., n], g, g0, dg, ndg, alfa = ConstantArray[0., n], d, c, m0, dm,
    lam, un, gn, it = 0, nav = 0, conv = False, hist = {}, beta},
   g = G[u]; nav++;
   g0 = If[g == 0, 1., Abs[g]];
   While[it < OptionValue["MaxIter"],
    it++;
    dg = gradiente[G, u, h]; nav += 2 n;
    ndg = Norm[dg]; alfa = -dg/ndg;
    AppendTo[hist, {it, u, g, Norm[u]}];
    If[Abs[g]/g0 < OptionValue["TolG"] && Norm[u - (alfa . u) alfa] < OptionValue["TolU"], conv = True; Break[]];
    d = ((dg . u - g)/ndg^2) dg - u;                   (* do ponto atual ao ponto HL-RF *)
    lam = 1.;
    If[TrueQ[OptionValue["Armijo"]],
     (* merito m(u) = |u|^2/2 + c|G(u)|, com c > |u|/|grad G| para d ser direcao de descida *)
     c = (2 Norm[u] + 1.)/ndg;
     m0 = u . u/2 + c Abs[g];
     dm = (u + c Sign[g] dg) . d;
     While[True,
      un = u + lam d; gn = G[un]; nav++;
      If[un . un/2 + c Abs[gn] <= m0 + 10.^-4 lam dm || lam < 2.^-10, Break[]];
      lam = lam/2],
     un = u + d; gn = G[un]; nav++];
    u = un; g = gn];
   beta = alfa . u;                                    (* com sinal: negativo se a origem esta na falha *)
   <|"beta" -> beta, "u" -> u, "alfa" -> alfa, "pf" -> CDF[normal01, -beta], "iteracoes" -> it, "avaliacoes" -> nav,
    "convergiu" -> conv, "historico" -> hist|>];


(* ::Section:: *)
(*3. Monte Carlo e amostragem por importancia*)


(* Gv: G vetorizada, recebe as amostras (uma por linha) e devolve a lista dos valores de G *)
monteCarlo[Gv_, n_Integer, nAm_Integer] := Module[{q, pf},
   q = N[Boole[Thread[Gv[RandomVariate[normal01, {nAm, n}]] <= 0]]];
   pf = Mean[q];
   <|"pf" -> pf, "cov" -> If[pf > 0, Sqrt[(1 - pf)/(nAm pf)], Infinity],
    "beta" -> If[0 < pf < 1, -InverseCDF[normal01, pf], Indeterminate], "q" -> q|>];

(* densidade de amostragem h = N(ustar, I): q = I(G <= 0) phi(v)/phi(v - ustar), com ln w = -v.ustar + |ustar|^2/2 *)
importancia[Gv_, ustar_?VectorQ, nAm_Integer] := Module[{v, q, pf},
   v = RandomVariate[normal01, {nAm, Length[ustar]}] + ConstantArray[ustar, nAm];
   q = N[Boole[Thread[Gv[v] <= 0]]] Exp[-v . ustar + ustar . ustar/2];
   pf = Mean[q];
   <|"pf" -> pf, "cov" -> If[pf > 0, StandardDeviation[q]/(Sqrt[nAm] pf), Infinity],
    "beta" -> If[0 < pf < 1, -InverseCDF[normal01, pf], Indeterminate], "q" -> q|>];

mediaAcumulada[q_List] := Accumulate[q]/Range[Length[q]];
(* pontos espacados em escala log para os graficos de convergencia *)
amostraLog[lista_List] := With[{k = Union[Round[10^Range[0., Log10[Length[lista]], 0.02]]]}, Transpose[{k, lista[[k]]}]];


(* ::Section:: *)
(*4. Problema 1: g = X1 X2 - X3 (o do codigo original)*)


mod1 = modelo[{NormalDistribution[40, 5], NormalDistribution[50, 2.5], NormalDistribution[1000, 200]}];
g1X[{x1_, x2_, x3_}] := x1 x2 - x3;            (* vale para um ponto ou para as colunas de muitas amostras *)
g1U[u_?VectorQ] := g1X[xDeU[mod1, u]];
g1Us[us_?MatrixQ] := g1X[Transpose[xDeUs[mod1, us]]];

(* referencia: P(X3 >= X1 X2), integrando nas duas primeiras variaveis *)
pfRef1 = NIntegrate[PDF[normal01, u1] PDF[normal01, u2] CDF[normal01, (1000 - (40 + 5 u1) (50 + 2.5 u2))/200],
   {u1, -9, 9}, {u2, -9, 9}, PrecisionGoal -> 8];
Print["referencia: pf = ", pfRef1, ", beta = ", -InverseCDF[normal01, pfRef1]];   (* 1.17688*10^-3 e 3.0415 *)

f1 = form[g1U, 3];
Print["FORM: beta = ", f1["beta"], ", pf = ", f1["pf"], ", u* = ", f1["u"], ", x* = ", xDeU[mod1, f1["u"]],
  ", alfa^2 = ", f1["alfa"]^2, "; ", f1["iteracoes"], " iteracoes, ", f1["avaliacoes"], " avaliacoes de g"];
(* beta = 3.04907, pf = 1.1477*10^-3 (2,5 % abaixo da referencia), u* = {-2.28993, -0.676678, 1.8961},
   x* = {28.55, 48.31, 1379.2}, alfa^2 = {0.564, 0.049, 0.387}; 6 iteracoes, 43 avaliacoes *)

mc1 = monteCarlo[g1Us, 3, 10^6];
is1 = importancia[g1Us, f1["u"], 10^3];
Print["Monte Carlo, N = 10^6: pf = ", mc1["pf"], " (CoV ", mc1["cov"], ")"];         (* CoV ~ 0,03 *)
Print["importancia, N = 10^3: pf = ", is1["pf"], " (CoV ", is1["cov"], ")"];         (* CoV ~ 0,06 *)
(* Monte Carlo com N = 1000 (o do codigo original): CoV ~ 0,92 e, em 31 % das rodadas, nenhuma falha (pf = 0) *)

qmc = monteCarlo[g1Us, 3, 10^5]["q"];
qis = importancia[g1Us, f1["u"], 10^5]["q"];
ListLogLinearPlot[{amostraLog[mediaAcumulada[qmc]], amostraLog[mediaAcumulada[qis]]}, Joined -> True,
 PlotRange -> {0, 0.003}, GridLines -> {None, {pfRef1, f1["pf"]}},
 PlotLegends -> {"Monte Carlo", "importancia (ponto de projeto)"}, AxesLabel -> {"N", "Pf acumulada"}]


(* ::Section:: *)
(*5. Problema 2 (proposto): talude infinito sob chuva*)


(* Plano de ruptura paralelo a superficie, a profundidade zf; Bishop com chi = 1 (phi^b = phi'), como no
   SlopeRandomField. h = carga de pressao no plano de ruptura (m): h < 0 succao, h > 0 pressao positiva.
   Parametros do talude de Vargas Ceron et al. (2025): beta = 45 graus, c' LN(10 kPa; 30 %), phi' LN(30 graus; 10 %).
   A chuva entra pela media mu_h de h: mais chuva, maior mu_h. *)
gama = 20.; gamaW = 9.81; betaS = 45. Degree; zf = 2.; sigH = 0.5;
sigN = gama zf Cos[betaS]^2; tau = gama zf Sin[betaS] Cos[betaS];    (* 20 kPa e 20 kPa *)
chi[h_] := 1;
(* alternativa: chi = Se(h) de van Genuchten, por exemplo com alfaVG = 1/m e nVG = 1.6:
   chi[h_] := UnitStep[h] + (1 - UnitStep[h]) (1 + (alfaVG Abs[h])^nVG)^(-(1 - 1/nVG)) *)
fsTalude[{c_, phi_, h_}] := (c + (sigN - gamaW h chi[h]) Tan[phi Degree])/tau;
g2X[x_] := fsTalude[x] - 1;

mod2[muh_, rho_: -0.5] := modelo[{lnDist[10., 0.3], lnDist[30., 0.1], NormalDistribution[muh, sigH]},
   {{1, rho, 0}, {rho, 1, 0}, {0, 0, 1}}];

(* referencia: dado phi' e h, c' e lognormal (condicionada a phi'), e a falha e c' < c'_crit(phi', h) *)
pfRef2[muh_, rho_: -0.5] := Module[{dc = lnDist[10., 0.3], dp = lnDist[30., 0.1], rz, cc},
   rz = rhoZ[dc, dp, rho];
   cc[z2_, h_] := tau - (sigN - gamaW h chi[h]) Tan[xDeZ[dp, z2] Degree];
   NIntegrate[PDF[normal01, z2] PDF[NormalDistribution[muh, sigH], h] *
     If[cc[z2, h] > 0, CDF[normal01, (Log[cc[z2, h]] - dc[[1]] - dc[[2]] rz z2)/(dc[[2]] Sqrt[1 - rz^2])], 0],
    {z2, -8, 8}, {h, muh - 8 sigH, muh + 8 sigH}, PrecisionGoal -> 8,
    Method -> {Automatic, "SymbolicProcessing" -> 0}]];

(* um cenario: referencia, FORM, Monte Carlo e importancia *)
cenario[muh_, rho_: -0.5] := Module[{m = mod2[muh, rho], G, Gv, f},
   G[u_?VectorQ] := g2X[xDeU[m, u]];
   Gv[us_?MatrixQ] := g2X[Transpose[xDeUs[m, us]]];
   f = form[G, 3];
   <|"muh" -> muh, "ref" -> pfRef2[muh, rho], "form" -> f, "x*" -> xDeU[m, f["u"]],
    "mc" -> monteCarlo[Gv, 3, 10^5], "is" -> importancia[Gv, f["u"], 10^3]|>];

(* cenario A (antes da chuva, mu_h = -2 m): ref 1.535*10^-4; FORM beta 3.5414, pf 1.99*10^-4; alfa^2 = {0.13, 0.18, 0.69}
   cenario B (depois da chuva, mu_h = -0.5 m): ref 0.1215 (o P_f do Ceron e 11,5 %); FORM beta 1.0781, pf 0.1405 *)
cA = cenario[-2.];
cB = cenario[-0.5];
Table[Print[c["muh"], ": ref ", c["ref"], "; FORM beta ", c["form"]["beta"], " pf ", c["form"]["pf"], " alfa^2 ",
   c["form"]["alfa"]^2, " x* ", c["x*"], "; MC ", c["mc"]["pf"], " (CoV ", c["mc"]["cov"], "); IS ", c["is"]["pf"],
   " (CoV ", c["is"]["cov"], ")"], {c, {cA, cB}}];

(* curva de fragilidade P_f(mu_h): referencia (rho = -0,5 e 0), FORM e importancia com 1000 amostras *)
frag = Table[Module[{m = mod2[mu], G, Gv, f, is},
    G[u_?VectorQ] := g2X[xDeU[m, u]];
    Gv[us_?MatrixQ] := g2X[Transpose[xDeUs[m, us]]];
    f = form[G, 3]; is = importancia[Gv, f["u"], 1000];
    {mu, pfRef2[mu], pfRef2[mu, 0], f["pf"], is["pf"], is["cov"], f["alfa"]^2}], {mu, Range[-4., 0.5, 0.25]}];
ListLogPlot[{frag[[All, {1, 2}]], frag[[All, {1, 3}]], frag[[All, {1, 4}]], frag[[All, {1, 5}]]},
 Joined -> {True, True, True, False}, PlotMarkers -> {None, None, None, Automatic},
 PlotLegends -> {"referencia, rho = -0,5", "referencia, rho = 0", "FORM", "importancia, N = 1000"},
 AxesLabel -> {"mu_h (m)", "Pf"}, PlotRange -> All]
ListLinePlot[Transpose[Table[{r[[1]], #} & /@ r[[7]], {r, frag}]], PlotLegends -> {"c'", "phi'", "h"},
 AxesLabel -> {"mu_h (m)", "alfa^2"}]


(* ::Section:: *)
(*6. As funcoes do codigo original, corrigidas*)


(* ::Text:: *)
(*Os mesmos nomes e argumentos; as mudancas estao marcadas. data = {{caso, media, dp}, ...}, caso 1 normal, 2 lognormal;*)
(*a normal equivalente de Rackwitz-Fiessler da lognormal depende do ponto x (e deve ser refeita a cada iteracao).*)


ComputeXZ[data_, x_] := Module[{caso, mu, sig, p1, p2},
   {caso, mu, sig} = data;
   Switch[caso,
    1, {sig, mu},
    2, p2 = Sqrt[Log[1 + (sig/mu)^2]]; p1 = Log[mu] - p2^2/2; {x p2, x (1 - Log[x] + p1)}] // N];
ComputeJXZ[data_, x_] := Module[{sz = Length[data], sn},
   sn = Table[ComputeXZ[data[[i]], x[[i]]], {i, sz}];               (* {dp, media} equivalentes *)
   {DiagonalMatrix[sn[[All, 1]]], DiagonalMatrix[1/sn[[All, 1]]], sn[[All, 2]]}];   (* {Jxz, Jzx, Mneq} *)
(* Y independentes, Z correlacionadas: z = L.y, y = L^-1.z, RX = L.Transpose[L] *)
ComputeJYZ[RX_] := Module[{L},
   L = Transpose[CholeskyDecomposition[N[RX]]];                     (* vale tambem para uma variavel: {{1.}} *)
   If[Max[Abs[L . Transpose[L] - RX]] > 10^-8, Print["Atencao: a decomposicao de Cholesky falhou!"]];   (* era Det *)
   {Inverse[L], L}];                                                 (* {Jyz, Jzy} *)
FromXtoY[data_, x_, RX_] := Module[{Jyz, Jzy, Jxz, Jzx, Mneq, Jyx, Jxy},
   {Jyz, Jzy} = ComputeJYZ[RX];
   {Jxz, Jzx, Mneq} = ComputeJXZ[data, x];
   Jyx = Jyz . Jzx;                                                  (* era Jzy.Jzx *)
   Jxy = Jxz . Jzy;                                                  (* era Jxz.Jyz *)
   {Jyx . (x - Mneq), Jxy, Jyx, Jxz, Mneq}];
FromYtoZ[y_, RX_] := Module[{Jyz, Jzy}, {Jyz, Jzy} = ComputeJYZ[RX]; Jzy . y];   (* nomes na ordem certa *)
FromZtoX[z_, Jxz_, Mneq_] := Jxz . z + Mneq;                         (* era Jxz.(z + Mneq) *)

(* teste de ida e volta com correlacao (com RX = I, o caso do codigo original, o erro nao aparece) *)
datarv = {{1, 40, 5}, {1, 50, 2.5}, {1, 1000, 200}};
RXc = {{1, 0.5, 0}, {0.5, 1, 0}, {0, 0, 1}};
x0 = {35., 52., 1100.};
{y0, Jxy0, Jyx0, Jxz0, Mneq0} = FromXtoY[datarv, x0, RXc];
Print["ida e volta x -> y -> z -> x: ", FromZtoX[FromYtoZ[y0, RXc], Jxz0, Mneq0], " (deve dar ", x0, ")"];


(* ::Section:: *)
(*7. Esboco: G de um modelo externo (o talude de elementos finitos)*)


(* ::Text:: *)
(*No talude com campos aleatorios, U sao os coeficientes xi da expansao de Karhunen-Loeve (ja normais padrao e*)
(*independentes) e G(u) = Gamma(u) - 1, com Gamma o fator do GIM. Um G externo teria a forma abaixo (nao executado:*)
(*o programa e o arquivo sao ilustrativos). Cada avaliacao e uma analise; o gradiente custa 2n analises.*)


(*
gTalude[u_?VectorQ] := Module[{arq = "xi.csv", saida},
   Export[arq, u];                                                   (* os xi da realizacao *)
   saida = RunProcess[{"./SlopeMonteCarlo", "xi", arq}, "StandardOutput"];
   ToExpression[Last[StringSplit[saida]]] - 1];                      (* Gamma - 1 *)
*)
