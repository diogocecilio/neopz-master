(* Substitutos das funcoes estatisticas que o Mathics nao tem (normal, lognormal, Cholesky, RandomVariate), so para o
   teste do confiabilidade.wl sem o Mathematica. NIntegrate e os graficos viram stubs: as referencias sao as do Python. *)
Unprotect[CDF, PDF, InverseCDF, Mean, StandardDeviation, RandomVariate, CholeskyDecomposition, NIntegrate, ListLogLinearPlot, ListLogPlot, ListLinePlot];
CDF[NormalDistribution[m_, s_], x_?NumericQ] := 0.5 Erfc[-(x - m)/(s Sqrt[2.])];
PDF[NormalDistribution[m_, s_], x_?NumericQ] := Exp[-((x - m)/s)^2/2]/(s Sqrt[2. Pi]);
InverseCDF[NormalDistribution[m_, s_], p_?NumericQ] := m + s Sqrt[2.] InverseErf[2 p - 1];
Mean[NormalDistribution[m_, s_]] := m;
StandardDeviation[NormalDistribution[m_, s_]] := s;
Mean[LogNormalDistribution[m_, s_]] := Exp[m + s^2/2];
StandardDeviation[LogNormalDistribution[m_, s_]] := Exp[m + s^2/2] Sqrt[Exp[s^2] - 1];
RandomVariate[NormalDistribution[0, 1], {n_, k_}] := Module[{a, b},
   a = RandomReal[{0, 1}, {n, k}]; b = RandomReal[{0, 1}, {n, k}];
   Sqrt[-2. Log[a]] Cos[2. Pi b]];
CholeskyDecomposition[m_?MatrixQ] := Module[{n = Length[m], U},
   U = Table[0., {n}, {n}];
   Do[U[[i, i]] = Sqrt[m[[i, i]] - Sum[U[[k, i]]^2, {k, 1, i - 1}]];
      Do[U[[i, j]] = (m[[i, j]] - Sum[U[[k, i]] U[[k, j]], {k, 1, i - 1}])/U[[i, i]], {j, i + 1, n}], {i, n}];
   U];
NIntegrate[___] := 0.00117688;
ListLogLinearPlot[___] := Null; ListLogPlot[___] := Null; ListLinePlot[___] := Null;
