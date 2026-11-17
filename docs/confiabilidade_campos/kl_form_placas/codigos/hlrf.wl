HLRF[medias_, Diag_, muvec_, tipos_, L_, gfunc_, gradfunc_, 
   maxiter_ : 5] := 
  Block[{gx, gradx, i, x, y, z, xk = medias, grady, checkconv, alpha, 
    beta, Jxy, Jyx, n = Length[tipos], Jyz = Inverse[L], Jzy = L},
   x = xk;
   z = Table[If[tipos[[j]] == "NORMAL", x[[j]], Log[x[[j]]]], {j, n}];
   y = Jyz . Inverse[Diag] . (z - muvec);
   For[i = 1, i <= maxiter, i++, 
    Jxy = DiagonalMatrix[
       Table[If[tipos[[j]] == "NORMAL", 1, x[[j]]], {j, n}]] . Diag . 
      Jzy;
    Jyx = Inverse[Jxy];
    gx = gfunc[x];
    gradx = gradfunc[x];
    grady = Transpose[Jxy] . gradx;
    y = ((grady . y - gx)/(grady . grady)) grady;
    z = muvec + Diag . L . y;
    x = Table[If[tipos[[j]] == "NORMAL", z[[j]], Exp[z[[j]]]], {j, n}];
    gx = gfunc[x];
    checkconv = Abs[gx];
    Jxy = 
     DiagonalMatrix[
       Table[If[tipos[[j]] == "NORMAL", 1, x[[j]]], {j, n}]] . Diag . 
      Jzy;
    gradx = gradfunc[x];
    grady = Transpose[Jxy] . gradx;
    alpha = (grady/Sqrt[grady . grady])^2;
    beta = Norm[y];
    Print[" | iter = ", i, " | residual = ", checkconv, " | G(x) = ", 
     N[gx], " | beta = ", beta];
    ];
   <|"x" -> x, "y" -> y, 
    "beta" -> beta, "alpha2Y" -> alpha, "residual" -> checkconv|>
   ];
