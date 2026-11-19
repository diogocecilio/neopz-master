(* nystrom_autopares.wl: os M maiores autopares de W^(1/2) C W^(1/2) da base de Nystrom do talude do Cho (2010), para
   o nystrom_base.py (opcao autopares=), com o Eigensystem do Mathematica (MKL, em paralelo): o eigh do Python sem
   threads leva horas com 15 000 nos.
   Uso: wolframscript -file nystrom_autopares.wl <nos> <M> <lh> <lv> <w> <prefixo>
     <nos>      os nos x y (um por linha), do nystrom_base.py nos=
     <w>        o peso de cada no (dx dy: a grade e uniforme)
   Grava <prefixo>.lam (os M autovalores, decrescentes) e <prefixo>.f (f_k(x_j) = v_k(x_j) / Sqrt[w], n x M, linha por
   no), em Real64. *)
{arqNos, mTxt, lhTxt, lvTxt, wTxt, pref} = Take[$ScriptCommandLine, -6];
{M, lh, lv, w} = {ToExpression[mTxt], N[ToExpression[lhTxt]], N[ToExpression[lvTxt]], N[ToExpression[wTxt]]};
pts = N[Import[arqNos, "Table"]];
n = Length[pts];
Print[n, " nos, M = ", M, ", lh = ", lh, ", lv = ", lv, ", w = ", w];
t = AbsoluteTiming[
   K = w Exp[-Abs[Outer[Subtract, pts[[All, 1]], pts[[All, 1]]]]/lh - Abs[Outer[Subtract, pts[[All, 2]], pts[[All, 2]]]]/lv];
   {lam, v} = Eigensystem[K, M];
   ClearAll[K]][[1]];
ordem = Ordering[lam, All, Greater];
lam = lam[[ordem]]; v = v[[ordem]];
Print["autopares em ", Round[t], " s; lambda_1 = ", lam[[1]], ", lambda_M = ", lam[[-1]], "; variancia retida = ",
  Total[lam]/(n w)];
f = Transpose[v]/Sqrt[w];                             (* n x M: f_k(x_j), linha por no *)
BinaryWrite[pref <> ".lam", lam, "Real64", ByteOrdering -> -1]; Close[pref <> ".lam"];
BinaryWrite[pref <> ".f", Flatten[f], "Real64", ByteOrdering -> -1]; Close[pref <> ".f"];
Print["gravados ", pref, ".lam e ", pref, ".f"];
