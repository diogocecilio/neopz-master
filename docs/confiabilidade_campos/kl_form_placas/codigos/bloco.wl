(* ===== cole DEPOIS do bloco Nystrom (precisa de mesh, coords, els, sz, w, KN, pars e das funcoes do FEM) ===== *)

(* 1. autovalores e autofuncoes nos nos *)
{lam, vec} = Eigensystem[KN];
ord = Ordering[lam, All, Greater]; lam = lam[[ord]]; vec = vec[[ord]];
phi = Table[vec[[i]]/Sqrt[w], {i, sz}];

(* 2. rigidez com E nodal variavel *)
CalcStiffE[elcoords_, eldisp_, Eel_, shapeU_, intrule_] := Block[{ek = 0., ef = 0., kp, fp},
   Do[Block[{pars = Append[pars, "young" -> shapeU[[ip, 1]] . Eel]},
      {kp, fp} = ContributeElastic[{shapeU[[ip, 1]], shapeU[[ip, 2]], elcoords, intrule[[ip, 4]], eldisp}]];
    ek += kp; ef += fp, {ip, Length[intrule]}];
   {ek, ef}];
AssembleE[mesh_, disp_, Enod_] := Block[{eltype, order, nvertex, intrule, shapeU, topol, xy, nU, d, ek, ef, bK, bF},
   xy = mesh["Coordinates"];
   {eltype, order, nvertex} = MeshInfo[mesh];
   intrule = IntegrationRule[eltype, pars["intdegree"]];
   shapeU = ShapeAtPoints[eltype, order, intrule];
   topol = mesh["MeshElements"][[1, 1]];
   nU = pars["dim"] Length[xy];
   {bK, bF} = Transpose@Table[d = ElementDofs[el];
      {ek, ef} = CalcStiffE[xy[[el]], disp[[d]], Enod[[el]], shapeU, intrule];
      {{d, d, ek}, {d, ef}}, {el, topol}];
   {GlobalMatrix[bK, {nU, nU}], GlobalVector[bF, nU]}];

(* 3. K0, K^(i), carga e contorno *)
nU = pars["dim"] sz;
zero = ConstantArray[0., nU];
K0 = First@AssembleE[mesh, zero, ConstantArray[1., sz]];
Kk = Table[First@AssembleE[mesh, zero, Sqrt[lam[[i]]] phi[[i]]], {i, sz}];
Fext = TopLoad[mesh, 1.];                                  (* tracao unitaria no topo *)
base = Flatten[Position[coords, {_, y_} /; Abs[y + 0.5] < 10^-8]];
free = Complement[Range[nU], ElementDofs[base]];
nA = First@Flatten@Position[coords, {x_, y_} /; Norm[{x, y} - {0.5, 0.5}] < 10^-8];
posA = First@First@Position[free, 2 nA];                   (* u_y no canto A *)
K0f = Normal[K0][[free, free]];
Kkf = Normal[#][[free, free]] & /@ Kk;
Ff = Fext[[free]];

(* 4. Monte Carlo *)
SeedRandom[1234];
nMC = 5000;
xi = RandomVariate[NormalDistribution[], {nMC, sz}];
MCrun[sig_?NumericQ, M_Integer] := Module[{u},
   u = Table[LinearSolve[K0f + sig Total[xi[[s, 1 ;; M]] Kkf[[1 ;; M]]], Ff][[posA]], {s, nMC}];
   {Mean[u], StandardDeviation[u]}];

sigs = Range[0., 0.30, 0.02];
res = Table[{s, MCrun[s, 2][[2]], MCrun[s, 4][[2]], MCrun[s, sz][[2]]}, {s, sigs}];
TableForm[res, TableHeadings -> {None, {"sigma_E", "M=2", "M=4", "M=25"}}]

book = {{0.10, 0.21}, {0.12, 0.26}, {0.14, 0.32}, {0.16, 0.37}, {0.18, 0.43}, {0.20, 0.48},
   {0.22, 0.55}, {0.24, 0.62}, {0.26, 0.70}, {0.28, 0.84}, {0.30, 0.99}};
Show[ListLinePlot[{res[[All, {1, 2}]], res[[All, {1, 3}]], res[[All, {1, 4}]]}, PlotMarkers -> Automatic,
   PlotLegends -> {"Nystrom M=2", "Nystrom M=4", "Nystrom M=25"}],
  ListPlot[book, PlotStyle -> Black, PlotMarkers -> {"*", 18}, PlotLegends -> {"Livro (Fig. 5.46)"}],
  Frame -> True, FrameLabel -> {"desvio padrao de E", "desvio de u_y(A)"}, PlotRange -> {{0, 0.31}, {0, 1.05}}]
