Needs["NDSolve`FEM`"];
pars = <|"young" -> 1., "nu" -> 0., "dim" -> 2, "intdegree" -> 2, "q" -> -1.|>;

IntegrationRule[eltype_, degree_] := Block[{pts, ws, rule},
   pts = ElementIntegrationPoints[eltype, degree];
   ws = ElementIntegrationWeights[eltype, degree];
   rule = Table[Join[PadRight[pts[[i]], 3], {ws[[i]]}], {i, Length[ws]}];
   rule];
ComputeShape[eltype_, order_] := Block[{psis, dpsis, xi, eta, zeta, vars},
   vars = Which[eltype === LineElement, {xi},
     eltype === QuadElement || eltype === TriangleElement, {xi, eta}, True, {xi, eta, zeta}];
   psis = ElementShapeFunction[eltype, order] @@ vars;
   dpsis = Table[D[psis, v], {v, vars}];
   {psis, dpsis}];
ShapeAtPoints[eltype_, order_, intrule_] := Block[{psis, dpsis, xi, eta, zeta, table},
   {psis, dpsis} = ComputeShape[eltype, order];
   table = Table[{xi, eta, zeta} = intrule[[i, 1 ;; 3]]; N[{psis, dpsis}], {i, Length[intrule]}];
   table];
ComputeBN[psis_, GradPhi_] := Block[{dim = pars["dim"], g = GradPhi, z, rows, Bu, Nu},
   z = ConstantArray[0., Length[psis]];
   rows = If[dim == 2, {{g[[1]], z}, {g[[2]], g[[1]]}, {z, z}, {z, g[[2]]}, {z, z}, {z, z}},
     {{g[[1]], z, z}, {g[[2]], g[[1]], z}, {g[[3]], z, g[[1]]}, {z, g[[2]], z}, {z, g[[3]], g[[2]]}, {z, z, g[[3]]}}];
   Bu = Flatten[Transpose[#]] & /@ rows;
   Nu = Table[Flatten[Transpose[Table[If[k == j, psis, z], {k, dim}]]], {j, dim}];
   {Bu, Nu}];
ContributeElastic[data_] := Block[{young = pars["young"], nu = pars["nu"], psis, GradPsi,
    elcoords, weight, eldisp, Jac, Bu, Nu, dvol, lam, mu, Cep, stress, K, Fint},
   {psis, GradPsi, elcoords, weight, eldisp} = data;
   Jac = GradPsi . elcoords;
   {Bu, Nu} = ComputeBN[psis, Inverse[Jac] . GradPsi];
   dvol = weight Det[Jac];
   lam = young nu/((1 + nu) (1 - 2 nu));
   mu = young/(2 (1 + nu));
   Cep = {{lam + 2 mu, 0, 0, lam, 0, lam}, {0, mu, 0, 0, 0, 0}, {0, 0, mu, 0, 0, 0},
     {lam, 0, 0, lam + 2 mu, 0, lam}, {0, 0, 0, 0, mu, 0}, {lam, 0, 0, lam, 0, lam + 2 mu}};
   stress = Cep . (Bu . eldisp);
   K = Transpose[Bu] . Cep . Bu dvol;
   Fint = Transpose[Bu] . stress dvol;
   {K, Fint}];

(* ---------- funcoes que faltavam ---------- *)
MeshInfo[mesh_] := Block[{elems = mesh["MeshElements"][[1]]},
   {Head[elems], mesh["MeshOrder"], Length[elems[[1, 1]]]}];
ElementDofs[el_] := Flatten[Table[pars["dim"] (n - 1) + k, {n, el}, {k, pars["dim"]}]];
CalcStiff[elcoords_, eldisp_, shapeU_, intrule_] := Block[{ek, ef, kp, fp},
   ek = 0.; ef = 0.;
   Do[{kp, fp} = ContributeElastic[{shapeU[[ip, 1]], shapeU[[ip, 2]], elcoords,
        intrule[[ip, 4]], eldisp}];
    ek += kp; ef += fp, {ip, Length[intrule]}];
   {ek, ef}];
GlobalMatrix[blocks_, dims_] := SparseArray[Normal@Merge[Flatten[
      Thread[Flatten[Outer[List, #[[1]], #[[2]]], 1] -> Flatten[#[[3]]]] & /@ blocks], Total], dims];
GlobalVector[blocks_, n_] := Block[{f = ConstantArray[0., n]},
   Do[f[[b[[1]]]] += b[[2]], {b, blocks}]; f];

Assemble[mesh_, disp_] := Block[{i, el, eltype, order, nvertex, intrule, shapeU, topol,
    coords, nU, ek, ef, d, bK, bF, K, Fint},
   coords = mesh["Coordinates"];
   {eltype, order, nvertex} = MeshInfo[mesh];
   intrule = IntegrationRule[eltype, pars["intdegree"]];
   shapeU = ShapeAtPoints[eltype, order, intrule];
   topol = mesh["MeshElements"][[1, 1]];
   nU = pars["dim"] Length[coords];
   bK = Table[0, {Length[topol]}]; bF = bK;
   For[i = 1, i <= Length[topol], i++, el = topol[[i]];
    d = ElementDofs[el];
    {ek, ef} = CalcStiff[coords[[el]], disp[[d]], shapeU, intrule];
    bK[[i]] = {d, d, ek}; bF[[i]] = {d, ef};];
   K = GlobalMatrix[bK, {nU, nU}];
   Fint = GlobalVector[bF, nU];
   {K, Fint}];

(* carga distribuida q (direcao y) nas arestas da face superior y = 0.5 *)
TopLoad[mesh_, q_] := Block[{coords = mesh["Coordinates"], edges, top, F, L},
   F = ConstantArray[0., pars["dim"] Length[coords]];
   edges = Join @@ (#[[1]] & /@ mesh["BoundaryElements"]);
   top = Select[edges, Abs[coords[[#[[1]], 2]] - 0.5] < 10^-8 && Abs[coords[[#[[2]], 2]] - 0.5] < 10^-8 &];
   Do[L = Norm[coords[[e[[2]]]] - coords[[e[[1]]]]];
    F[[2 e]] += q L/2, {e, top}];
   F];

(* ---------- problema ---------- *)
mesh = ToElementMesh[FullRegion[2], {{-0.5, 0.5}, {-0.5, 0.5}}, "MeshOrder" -> 1, MaxCellMeasure -> 0.1];
coords = mesh["Coordinates"];
nU = pars["dim"] Length[coords];
disp = ConstantArray[0., nU];
{K, Fint} = Assemble[mesh, disp];
Fext = TopLoad[mesh, pars["q"]];
base = Flatten[Position[coords, {_, y_} /; Abs[y + 0.5] < 10^-8]];
fixed = ElementDofs[base];
free = Complement[Range[nU], fixed];
disp[[free]] = LinearSolve[K[[free, free]], (Fext - Fint)[[free]]];
ux = disp[[1 ;; ;; 2]]; uy = disp[[2 ;; ;; 2]];
Print[Head[mesh["MeshElements"][[1]]]];Print["Total carga: ", Total[Fext]];
Print["max |ux| = ", Max[Abs[ux]], "   min uy = ", Min[uy]];
Print["erro vs exato uy=-(y+0.5): ", Max[Abs[uy + (coords[[All, 2]] + 0.5)]]];

uxI = ElementMeshInterpolation[{mesh}, ux];
uyI = ElementMeshInterpolation[{mesh}, uy];
pUy = ContourPlot[uyI[x, y], {x, y} \[Element] mesh, ColorFunction -> "TemperatureMap",
   Contours -> 15, PlotLegends -> Automatic, PlotLabel -> "uy", AspectRatio -> Automatic];
pUx = ContourPlot[uxI[x, y], {x, y} \[Element] mesh, ColorFunction -> "TemperatureMap",
   Contours -> 15, PlotLegends -> Automatic, PlotLabel -> "ux", AspectRatio -> Automatic];
scale = 0.2;
meshDef = ToElementMesh["Coordinates" -> coords + scale Transpose[{ux, uy}],
   "MeshElements" -> mesh["MeshElements"]];
pDef = Show[mesh["Wireframe"["MeshElementStyle" -> Directive[FaceForm[], EdgeForm[Gray]]]],
   meshDef["Wireframe"["MeshElementStyle" -> Directive[FaceForm[], EdgeForm[{Red, Thick}]]]],
   PlotLabel -> "Deformada (escala 0.2)", Frame -> True];
Export["campos.png", GraphicsRow[{pUx, pUy, pDef}, ImageSize -> 1200]];
