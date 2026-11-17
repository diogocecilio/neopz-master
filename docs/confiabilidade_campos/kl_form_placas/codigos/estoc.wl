Get["base.wl"];
pars["q"] = 1.;   (* tracao, como no livro (plate stretching) *)

TopLoad[mesh_, q_] := Block[{coords = mesh["Coordinates"], edges, top, F, L},
   F = ConstantArray[0., pars["dim"] Length[coords]];
   edges = Join @@ (#[[1]] & /@ mesh["BoundaryElements"]);
   top = Select[edges, Abs[coords[[#[[1]], 2]] - 0.5] < 10^-8 && Abs[coords[[#[[2]], 2]] - 0.5] < 10^-8 &];
   Do[L = Norm[coords[[e[[2]]]] - coords[[e[[1]]]]]; F[[2 e]] += q L/2, {e, top}];
   F];

(* E variavel: interpolado dos nos para o ponto de integracao; reaproveita ContributeElastic *)
CalcStiffE[elcoords_, eldisp_, Eel_, shapeU_, intrule_] := Block[{ek = 0., ef = 0., kp, fp},
   Do[Block[{pars = Append[pars, "young" -> shapeU[[ip, 1]] . Eel]},
      {kp, fp} = ContributeElastic[{shapeU[[ip, 1]], shapeU[[ip, 2]], elcoords, intrule[[ip, 4]], eldisp}]];
    ek += kp; ef += fp, {ip, Length[intrule]}];
   {ek, ef}];
AssembleE[mesh_, disp_, Enod_] := Block[{eltype, order, nvertex, intrule, shapeU, topol, coords, nU, d, ek, ef, bK, bF},
   coords = mesh["Coordinates"];
   {eltype, order, nvertex} = MeshInfo[mesh];
   intrule = IntegrationRule[eltype, pars["intdegree"]];
   shapeU = ShapeAtPoints[eltype, order, intrule];
   topol = mesh["MeshElements"][[1, 1]];
   nU = pars["dim"] Length[coords];
   {bK, bF} = Transpose@Table[d = ElementDofs[el];
      {ek, ef} = CalcStiffE[coords[[el]], disp[[d]], Enod[[el]], shapeU, intrule];
      {{d, d, ek}, {d, ef}}, {el, topol}];
   {GlobalMatrix[bK, {nU, nU}], GlobalVector[bF, nU]}];

mesh = ToElementMesh[FullRegion[2], {{-0.5, 0.5}, {-0.5, 0.5}}, "MeshOrder" -> 1, MaxCellMeasure -> 0.1];

(* ---------- Nystrom (codigo do usuario) ---------- *)
coords = mesh["Coordinates"];
els = mesh["MeshElements"][[1, 1]];
sz = Length[coords];
LX = 1.; LY = 1.;
kernel[x1_, x2_] := Exp[-Abs[x1[[1]] - x2[[1]]]/LX - Abs[x1[[2]] - x2[[2]]]/LY];
CMAT = Table[Table[kernel[coords[[i]], coords[[j]]], {i, 1, sz}], {j, 1, sz}];
w = ConstantArray[0., sz];
Do[w[[t]] += Area[Polygon[coords[[t]]]]/Length[t], {t, els}];
Print["Total[w] = ", Total[w]];
sw = DiagonalMatrix[Sqrt[w]];
KN = sw . CMAT . sw;
(* ---------- autovalores e autofuncoes nos nos ---------- *)
{lam, vec} = Eigensystem[KN];
ord = Ordering[lam, All, Greater]; lam = lam[[ord]]; vec = vec[[ord]];
phi = Table[vec[[k]]/Sqrt[w], {k, sz}];          (* f_k(x_n) = v_k/sqrt(w_n) *)
Print["ortonormalidade max erro: ", Max@Abs[Table[Total[w phi[[i]] phi[[j]]], {i, 6}, {j, 6}] - IdentityMatrix[6]]];

(* ---------- KL analitico (livro, eq. 5.56-5.59), a = 0.5, c = 1 ---------- *)
a = 0.5; c = 1.;
om1 = Table[If[OddQ[i], x /. FindRoot[c - x Tan[x a], {x, (i - 1) Pi/(2 a) + 0.3}],
     x /. FindRoot[x + c Tan[x a], {x, (i - 1) Pi/(2 a) + 0.3}]], {i, 8}];
lam1 = 2 c/(om1^2 + c^2);
lamAn = Take[Sort[Flatten[Outer[Times, lam1, lam1]], Greater], 8];
Print["lambda Nystrom : ", NumberForm[Take[lam, 8], 4]];
Print["lambda analitico: ", NumberForm[lamAn, 4]];
Print["soma lambda (deve ser |A|=1): ", Total[lam]];

(* ---------- matrizes K0 e K^(k) (eq. 5.72 do livro) ---------- *)
nU = pars["dim"] sz; zero = ConstantArray[0., nU];
K0 = First@AssembleE[mesh, zero, ConstantArray[1., sz]];
Kk = Table[First@AssembleE[mesh, zero, Sqrt[lam[[k]]] phi[[k]]], {k, sz}];
Fext = TopLoad[mesh, pars["q"]];
base = Flatten[Position[coords, {_, y_} /; Abs[y + 0.5] < 10^-8]];
free = Complement[Range[nU], ElementDofs[base]];
nA = First@Flatten@Position[coords, {x_, y_} /; Norm[{x, y} - {0.5, 0.5}] < 10^-8];
dofA = 2 nA;   (* uy no canto A *)
uDet = LinearSolve[Normal[K0][[free, free]], Fext[[free]]];
uAdet = (ReplacePart[zero, Thread[free -> uDet]])[[dofA]];
Print["u_y(A) deterministico = ", uAdet];

(* ---------- Monte Carlo ---------- *)
SeedRandom[1234];
nMC = 5000;
xi = RandomVariate[NormalDistribution[], {nMC, sz}];
K0f = Normal[K0][[free, free]]; Kkf = Normal[#][[free, free]] & /@ Kk; Ff = Fext[[free]];
posA = First@First@Position[free, dofA];
MCrun[sig_, M_] := Module[{u, Emin},
   u = Table[LinearSolve[K0f + sig Sum[xi[[s, k]] Kkf[[k]], {k, M}], Ff][[posA]], {s, nMC}];
   Emin = Table[Min[1 + sig Sum[xi[[s, k]] Sqrt[lam[[k]]] phi[[k]], {k, M}]], {s, nMC}];
   {Mean[u], StandardDeviation[u], Count[Emin, _?(# <= 0 &)]}];
sigs = Range[0., 0.30, 0.02];
book = {{0.10, 0.21}, {0.12, 0.26}, {0.14, 0.32}, {0.16, 0.37}, {0.18, 0.43}, {0.20, 0.48},
   {0.22, 0.55}, {0.24, 0.62}, {0.26, 0.70}, {0.28, 0.84}, {0.30, 0.99}};
res = Table[Join[{s}, MCrun[s, 2], MCrun[s, 4], MCrun[s, sz]], {s, sigs}];
Print["sig | M=2: media, desvio, nE<=0 | M=4: media, desvio, nE<=0"];
Do[Print[ToString[NumberForm[#, {4, 3}]] & /@ r], {r, res}];
Export["res.m", {res, book, lam, lamAn}];

(* variancia pontual do campo truncado (fracao capturada) *)
Print["var captada no ponto A, M=2: ", Sum[lam[[k]] phi[[k, nA]]^2, {k, 2}], "  M=4: ", Sum[lam[[k]] phi[[k, nA]]^2, {k, 4}],
  "  M=25: ", Sum[lam[[k]] phi[[k, nA]]^2, {k, sz}]];

p1 = ListLinePlot[{res[[All, {1, 3}]], res[[All, {1, 6}]], res[[All, {1, 9}]]},
   PlotMarkers -> Automatic, PlotLegends -> {"Nystrom MC, M=2", "Nystrom MC, M=4", "Nystrom MC, M=25 (todos)"}];
p2 = ListPlot[book, PlotStyle -> Black, PlotMarkers -> {"*", 18}, PlotLegends -> {"Livro: MCS 5000 (Fig. 5.46)"}];
g = Show[p1, p2, Frame -> True, FrameLabel -> {"desvio padrao de E", "desvio de u_y(A) / u_det"},
   PlotRange -> {{0, 0.31}, {0, 1.05}}, ImageSize -> 600, GridLines -> Automatic];
pf = Table[ListDensityPlot[Transpose[{coords[[All, 1]], coords[[All, 2]], phi[[k]]}],
    ColorFunction -> "TemperatureMap", PlotLabel -> Row[{"f", k, ", \[Lambda]=", NumberForm[lam[[k]], 3]}],
    PlotLegends -> Automatic, ImageSize -> 220], {k, 4}];
Export["comparacao.png", g];
Export["autofuncoes.png", GraphicsRow[pf, ImageSize -> 1000]];
