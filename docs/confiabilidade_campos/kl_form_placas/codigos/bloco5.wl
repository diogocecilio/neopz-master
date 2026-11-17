(* ===== FORM com X = valores nodais de E (campo discretizado ponto a ponto) ===== *)
(* K e linear em E nodal: K(E) = Sum_n E_n Kn_n, com Kn_n = rigidez de uma "funcao chapeu" no no n *)
Knod = Table[Normal[First@AssembleE[mesh, zero, UnitVector[sz, n]]][[free, free]], {n, sz}];
SolveGradE[Ev_] := Module[{LS, u, lamb},
   LS = LinearSolve[Total[Ev Knod]];
   u = LS[Ff]; lamb = LS[UnitVector[Length[Ff], posA]];
   {u[[posA]], Table[-lamb . (Knod[[n]] . u), {n, sz}]}];   (* d u_A / d E_n *)

delta = 1.4;
GXE[Ev_] := delta - First@SolveGradE[Ev];
GRADXE[Ev_] := -Last@SolveGradE[Ev];

tiposE = ConstantArray["NORMAL", sz];
mediasE = ConstantArray[1., sz];          (* <E> = 1 em cada no *)
desviosE = ConstantArray[sigE, sz];       (* sigma_E em cada no *)
corrE = CMAT;                             (* correlacao entre nos = kernel exponencial *)
LE = Transpose[CholeskyDecomposition[corrE]];
rE = Block[{Print}, HLRF[mediasE, DiagonalMatrix[desviosE], mediasE, tiposE, LE, GXE, GRADXE, 10]];
Print["X = E nodal, Cholesky: beta = ", rE["beta"], "  Pf = ", CDF[NormalDistribution[], -rE["beta"]]];
Print["xi do KL, Nystrom M=25: beta = ", FORMrun[delta, sigE, sz, Kkf]["beta"]];
Print["E no ponto de projeto (nos do topo): ", rE["x"][[Flatten@Position[coords, {_, y_} /; y > 0.49]]]];

(* versao LOGNORMAL do mesmo campo *)
p2[mean_, sdev_] := Sqrt[Log[1 + (sdev/mean)^2]];
p1[mean_, sigma_] := Log[mean] - 0.5 sigma^2;
tiposL = ConstantArray["LOGNORMAL", sz];
sigL = p2[#, sigE] & /@ mediasE; muL = MapThread[p1, {mediasE, sigL}];
rL = Block[{Print}, HLRF[mediasE, DiagonalMatrix[sigL], muL, tiposL, LE, GXE, GRADXE, 10]];
Print["campo LOGNORMAL: beta = ", rL["beta"], "  Pf = ", CDF[NormalDistribution[], -rL["beta"]], " residuo ", rL["residual"]];
