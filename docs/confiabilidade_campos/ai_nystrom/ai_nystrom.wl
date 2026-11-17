(* ai_nystrom.wl: amostragem por importancia (AI) com a KL de Nystrom no talude do Cho (2010).
   O FS do elementos finitos e trocado por um estado limite barato (a media de c numa banda de ruptura), para tudo
   rodar em segundos no Mathematica. A parte da AI e a mesma do Cho2010 no NeoPZ. *)

Needs["NDSolve`FEM`"];
SeedRandom[2010];
aqui = If[$InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];

(* ================= 1. a KL de Nystrom (como no nymstron.nb) ================= *)
reg = Polygon[{{0, 0}, {30, 0}, {30, 5}, {20, 5}, {10, 15}, {0, 15}}];   (* o talude do Cho *)
mesh = ToElementMesh[reg, "MeshOrder" -> 1, MaxCellMeasure -> 0.25];
coords = mesh["Coordinates"];
els = mesh["MeshElements"][[1, 1]];
sz = Length[coords];
LX = 20.; LY = 2.;                                         (* comprimentos de correlacao *)
x = coords[[All, 1]]; y = coords[[All, 2]];
CMAT = Exp[-Abs[Outer[Subtract, x, x]]/LX - Abs[Outer[Subtract, y, y]]/LY];
w = ConstantArray[0., sz];                                 (* pesos: area dos elementos dividida entre os nos *)
Do[w[[t]] += Area[Polygon[coords[[t]]]]/Length[t], {t, els}];
sq = Sqrt[w];
KN = CMAT*Outer[Times, sq, sq];
M = 50;                                                    (* termos da KL *)
{lam, vec} = Eigensystem[KN, M];
modos = Table[vec[[k]]/sq, {k, 1, M}];                    (* f_k nos nos *)
Print["nos: ", sz, "   area: ", Total[w], "   variancia retida com M = ", M, ": ", Total[lam]/Total[w]];

(* ================= 2. o campo de c' (lognormal) ================= *)
mc = 10.; covc = 0.3;                                      (* media e CoV de c' (kPa), como no Cho *)
sln = Sqrt[Log[1 + covc^2]]; mln = Log[mc] - sln^2/2;
(* o campo nos nos para os M numeros xi: Y = soma dos M termos sqrt(lam_k) f_k xi_k; c = exp(mln + sln Y) *)
campo[xi_] := Exp[mln + sln Sum[Sqrt[lam[[k]]] modos[[k]] xi[[k]], {k, 1, M}]]

(* ================= 3. a banda de ruptura e a falha ================= *)
xc = 22.; yc = 24.; R = Norm[{20., 5.} - {xc, yc}];       (* um circulo pelo pe do talude *)
banda = {};
For[j = 1, j <= sz, j++,
 If[Abs[Norm[coords[[j]] - {xc, yc}] - R] < 1.,
   AppendTo[banda, j];
   ];
 ]
wb = Table[w[[banda[[j]]]], {j, 1, Length[banda]}];
(* a parte do campo so nos nos da banda: Bb[[k]] = sqrt(lam_k) f_k na banda (para ir rapido) *)
Bb = Table[Sqrt[lam[[k]]] modos[[k, banda]], {k, 1, M}];
cmedio[xi_] := Total[wb Exp[mln + sln xi . Bb]]/Total[wb]   (* media de c na banda *)
(* a direcao da banda: quanto cada xi_k muda, em media, o log de c na banda *)
a = Table[Total[wb Bb[[k]]]/Total[wb], {k, 1, M}];
d = -a/Norm[a];                                            (* andar em d enfraquece a banda *)
beta0 = 2.5;
ccrit = Exp[mln - beta0 sln Norm[a]];                      (* falha: media de c na banda < ccrit *)
falhou[xi_] := cmedio[xi] < ccrit
Print["nos na banda: ", Length[banda], "   c critico: ", ccrit];

(* ================= 4. Monte Carlo direto ================= *)
nd = 100000;
faild = 0;
For[i = 1, i <= nd, i++,
 xi = RandomVariate[NormalDistribution[0, 1], M];
 If[falhou[xi], faild++];
 ]
pfD = N[faild/nd]; covD = Sqrt[(1 - pfD)/(nd pfD)];
Print["direto: ", nd, " amostras, ", faild, " falhas, Pf = ", pfD, ", CoV = ", covD];

(* ================= 5. tres deslocamentos ================= *)
(* (a) banda + entropia cruzada: piloto direto de 1000; elite = as 10% de menor media de c; s = media de xi.d nela *)
np = 1000;
piloto = Table[RandomVariate[NormalDistribution[0, 1], M], {np}];
gp = Table[cmedio[piloto[[i]]], {i, 1, np}];
ordem = Ordering[gp];
nel = Max[Count[gp, v_ /; v < ccrit], Round[0.1 np]];
elite = Table[piloto[[ordem[[i]]]], {i, 1, nel}];
s = Mean[Table[elite[[i]] . d, {i, 1, nel}]];
muBanda = s d;
(* (b) os 10 primeiros modos: a media dos xi da elite nesses termos *)
K = 10;
muModos = Join[Mean[Table[elite[[i, 1 ;; K]], {i, 1, nel}]], ConstantArray[0., M - K]];
(* (c) FORM (HL-RF) em g = log(media de c / ccrit), gradiente por diferencas centrais *)
g[xi_] := Log[cmedio[xi]/ccrit]
u = ConstantArray[0., M]; h = 0.1;
For[it = 1, it <= 10, it++,
 g0 = g[u];
 grad = Table[(g[u + h UnitVector[M, k]] - g[u - h UnitVector[M, k]])/(2 h), {k, 1, M}];
 u = (grad . u - g0)/(grad . grad) grad;
 ]
betaFORM = Norm[u]; muFORM = u;
Print["FORM: beta = ", betaFORM, ", Pf = Phi(-beta) = ", CDF[NormalDistribution[], -betaFORM]];
Print["|mu|^2: banda ", muBanda . muBanda, ", modos ", muModos . muModos, ", FORM ", muFORM . muFORM];

(* ================= 6. a amostragem por importancia ================= *)
ai[mu_, n_] := Module[{fail = 0, pesos = {}, xi, pf, cv},
  For[i = 1, i <= n, i++,
   xi = RandomVariate[NormalDistribution[0, 1], M] + mu;   (* sorteio de h = N(mu, I) *)
   If[falhou[xi],
     AppendTo[pesos, Exp[-xi . mu + mu . mu/2]];            (* peso f/h *)
     fail++;
     ];
   ];
  pf = Total[pesos]/n;
  cv = Sqrt[(Total[pesos^2]/n - pf^2)/(n - 1)]/pf;
  {pf, cv, fail, Total[pesos]^2/Total[pesos^2], Max[pesos]/Total[pesos]}]
n = 2000;
res = Table[ai[mu, n], {mu, {muBanda, muModos, muFORM}}];
nomes = {"banda", "10 modos", "FORM"};
Do[Print["AI ", nomes[[k]], ": Pf = ", res[[k, 1]], ", CoV = ", res[[k, 2]], ", falhas = ", res[[k, 3]],
   ", ESS_F = ", res[[k, 4]], ", maior parcela = ", res[[k, 5]],
   ", ganho = ", (covD^2 nd)/(res[[k, 2]]^2 n)], {k, 3}];

(* ================= 7. figuras ================= *)
contorno = {Black, Line[{{0, 0}, {30, 0}, {30, 5}, {20, 5}, {10, 15}, {0, 15}, {0, 0}}]};
mapa[v_, titulo_, faixa_] := Module[{f = ElementMeshInterpolation[{mesh}, v]},
  DensityPlot[f[xx, yy], {xx, yy} \[Element] reg, ColorFunction -> "TemperatureMap", PlotPoints -> 50,
   PlotRange -> {All, All, faixa}, ColorFunctionScaling -> True, AspectRatio -> Automatic, Frame -> True,
   PlotLabel -> titulo, Epilog -> {contorno, Gray, PointSize[0.006], Point[coords[[banda]]]}, ImageSize -> 330]];
figModos = GraphicsGrid[Partition[Table[mapa[modos[[k]], "modo " <> ToString[k], All], {k, 6}], 3]];
Export[FileNameJoin[{aqui, "figs", "modos.png"}], figModos, ImageResolution -> 150];
figCampos = GraphicsGrid[{{mapa[campo[RandomVariate[NormalDistribution[0, 1], M]], "c' sorteado de f", {3, 18}],
    mapa[campo[RandomVariate[NormalDistribution[0, 1], M] + muFORM], "c' sorteado de h (FORM)", {3, 18}]}}];
Export[FileNameJoin[{aqui, "figs", "campos.png"}], figCampos, ImageResolution -> 150];
tmu = muFORM . d;
figDist = Plot[{PDF[NormalDistribution[0, 1], t], PDF[NormalDistribution[tmu, 1], t]}, {t, -4, 6},
   PlotLegends -> {"f: original", "h: deslocada (FORM)"}, AxesLabel -> {"t = xi.d", "densidade"},
   Epilog -> {Dashed, Line[{{beta0, 0}, {beta0, 0.42}}]}, ImageSize -> 420];
Export[FileNameJoin[{aqui, "figs", "distribuicoes.png"}], figDist, ImageResolution -> 150];
Print["figuras em ", FileNameJoin[{aqui, "figs"}]];
