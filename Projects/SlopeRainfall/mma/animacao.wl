(* ::Package:: *)
(* Pós-processamento de um caso do SlopeRainfall (talude não saturado sob chuva, Rahardjo et al. 2007):
   FS ao longo do tempo (redução de resistência e aumento da gravidade) e animação das isolinhas de p_w com o FS.
   Mesma estrutura do pós-processamento do rebaixamento (SlopeDrawdown); os dados vêm dos arquivos CSV do caso.

   Uso: wolframscript -file animacao.wl <pasta do caso>       (ou abra no Mathematica e defina dir)
   Saídas na pasta do caso: fs_tempo.png, chuva.mp4, chuva.gif *)

dir = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]],
   If[ValueQ[dir], dir, Directory[]]];
SetDirectory[dir];

(*----------------------------- dados do caso -----------------------------*)
xy = N@Rest@Import["vertices.csv", "CSV"];
tri = Rest@Import["triangles.csv", "CSV"] + 1;
cantos = N@Rest@Import["outline.csv", "CSV"];
snapsW = N@Import["snaps_pw.csv", "CSV"];          (* cada linha: t (h), p_w nos vértices *)
fsTab = N@Rest@Import["fs.csv", "CSV"];            (* t (h), FS SRM, FS GI, ... *)
summ = Association[Rule @@@ (StringSplit[#, Whitespace, 2] & /@ ReadList["summary.txt", String])];
num[k_] := ToExpression[summ[k]];
Ir = num["I"]; tRain = num["train"];
tFail = If[num["t_fail"] > 0, 3600. num["t_fail"], None];     (* s *)
{L, H} = Max /@ Transpose[cantos];

(* snaps: {t (s), I (mm/h) no instante, p_w nos vértices} *)
snaps = {3600. #[[1]], If[0 < #[[1]] <= tRain + 10^-9, Ir, 0.], Rest[#]} & /@ snapsW;

(* pares {t (h), FS}, t = 0 em 0,01 h; FS <= 0: não calculado *)
fsT = Select[{Max[#[[1]], 0.01], #[[2]]} & /@ fsTab, #[[2]] > 0 &];
fsGT = Select[{Max[#[[1]], 0.01], #[[3]]} & /@ fsTab, #[[2]] > 0 &];
If[fsGT === {}, fsGT = fsT];

tlabel[t_] := If[t < 3600., ToString[NumberForm[t/3600., {3, 2}]] <> " h", ToString[Round[t/3600.]] <> " h"];

(*----------------------------- FS em log(t) -----------------------------*)
(*FS interpolado em log(t) para o marcador da animação (constante se houver um só ponto)*)
InterpFS[l_] :=
  If[Length[l] < 2, Function[t, l[[1, 2]]],
   With[{f = Interpolation[{Log[#[[1]]], #[[2]]} & /@ l, InterpolationOrder -> 1], mm = MinMax[l[[All, 1]]]},
    Function[t, f[Log[Clip[Max[t/3600., 0.01], mm]]]]]];
fsAt = InterpFS[fsT];
fsGAt = InterpFS[fsGT];

fsPlot =
 ListLogLinearPlot[{fsT, fsGT}, Joined -> True,
  PlotMarkers -> {Automatic, 8},
  PlotStyle -> {Darker[Blue], Darker[Green]},
  PlotLegends -> Placed[{"redução de resistência", "aumento da gravidade"}, Below],
  Frame -> True,
  FrameLabel -> {"t (h)  [t = 0 desenhado em 0,01 h]", "FS"},
  PlotRange -> {All, {Min[0.9, Min[Join[fsT, fsGT][[All, 2]]] - 0.05], 1.05 Max[Join[fsT, fsGT][[All, 2]]]}},
  GridLines -> {{1.}, {1.}}, GridLinesStyle -> Directive[Gray, Dashed],
  Epilog -> If[tFail === None, {}, {Red, Thick, Line[{{Log[tFail/3600.], 0.}, {Log[tFail/3600.], 100.}}],
     Text[Style["ruptura", Red], {Log[tFail/3600.], 1.}, {-1.1, -1}]}],
  PlotLabel -> "Chuva de " <> ToString[Ir] <> " mm/h por " <> ToString[tRain] <> " h: FS ao longo do tempo",
  ImageSize -> 500];

Export["fs_tempo.png", fsPlot];

(*-----------------------------isolinhas de p_w e FS: animação-----------------------------*)
(*p_w vive nos vértices: malha linear dos triângulos do NeoPZ*)
meshP = ToElementMesh["Coordinates" -> xy, "MeshElements" -> {TriangleElement[tri]}];

pTop = 50. Ceiling[Max[snaps[[All, 3]]]/50.];
pRange = {-80., pTop};

(* chuva: setas sobre a superfície do terreno *)
superficie[x_] := Piecewise[{{3 H/3, x <= cantos[[5, 1]]},
    {H - (H/3) (x - cantos[[5, 1]])/(cantos[[4, 1]] - cantos[[5, 1]]), x <= cantos[[4, 1]]}}, 2 H/3];

frame[i_] :=
  Module[{t, ir, pv, f, th, slope, rain},
   {t, ir, pv} = snaps[[i, 1 ;; 3]];
   f = ElementMeshInterpolation[{meshP}, pv];
   th = Max[t/3600., 0.01];
   (*chuva e contorno do talude*)
   rain = Graphics[{If[ir > 0, {Darker[Blue], Table[Arrow[{{x, superficie[x] + 0.08 H}, {x, superficie[x] + 0.4}}],
         {x, 0.02 L, 0.98 L, 0.96 L/27}]}, {}],
      Black, Thin, Line[Append[cantos, First[cantos]]]}];
   slope =
    Legended[
     Show[(*preenchimento a cada 10 kPa*)
      ContourPlot[f[x, y], Element[{x, y}, meshP],
       Contours -> Range[-80, pTop, 10], ContourStyle -> None,
       ColorFunctionScaling -> False,
       ColorFunction -> (ColorData[{"TemperatureMap", pRange}][#] &)],(*linhas rotuladas a cada 50 kPa*)
      ContourPlot[f[x, y], Element[{x, y}, meshP],
       Contours -> Range[0, pTop, 50], ContourShading -> None,
       ContourStyle -> Directive[Black, Thin],
       ContourLabels -> All],(*p_w=0*)
      ContourPlot[f[x, y], Element[{x, y}, meshP], Contours -> {0.},
       ContourShading -> None,
       ContourStyle -> Directive[Magenta, Thick]], rain,
      Frame -> True, FrameLabel -> {"x (m)", "y (m)"},
      AspectRatio -> Automatic, PlotRange -> {{0, L}, {0, 1.13 H}},
      PlotLabel -> "t = " <> tlabel[t] <> "    " <> If[ir > 0, "chuva " <> ToString[ir] <> " mm/h", "sem chuva"] <>
        "    (magenta: p_w = 0)",
      ImageSize -> 560],
     BarLegend[{"TemperatureMap", pRange}, LegendLabel -> "p_w (kPa)"]];
   Row[{slope,
     Show[fsPlot,
      Graphics[{Red, Dashed, Line[{{Log[th], 0.}, {Log[th], 100.}}],
        PointSize[0.025], Point[{Log[th], fsAt[t]}],
        Point[{Log[th], fsGAt[t]}]}],
      PlotLabel -> "FS = " <> ToString[NumberForm[fsAt[t], {3, 2}]] <>
        "    FSg = " <> ToString[NumberForm[fsGAt[t], {3, 2}]],
      ImageSize -> 420]}]];

frames = Table[Rasterize[frame[i], ImageResolution -> 96], {i, Length[snaps]}];
(*mesmo tamanho em todos os quadros*)
frames = ImageResize[#, ImageDimensions[First[frames]]] & /@ frames;

If[$Notebooks, ListAnimate[frames, AnimationRunning -> False]];

Export["chuva.mp4", frames, "FrameRate" -> 4];
Export["chuva.gif", frames, "DisplayDurations" -> 0.25, "AnimationRepetitions" -> Infinity];
