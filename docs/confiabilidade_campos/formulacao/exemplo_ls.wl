(* exemplo_ls.wl: o line sampling (amostragem por linhas) em exemplos com P_f exata, para a Secao "Line sampling" da
   formulacao. Rodar: wolframscript -file exemplo_ls.wl (grava as figuras em figs_ls/ e imprime os numeros da saida).
   Partes:
   1. a figura do metodo (2D, superficie de falha curva);
   2. superficie curva em 2D: P_f exata, line sampling, amostragem por importancia (AI) e Monte Carlo direto, 200
      repeticoes de cada, com o mesmo numero de analises;
   3. superficie plana: o line sampling tem variancia zero;
   4. uma direcao errada (girada) nao cria vies, so aumenta a variancia;
   5. muitas dimensoes (n = 2000, um "campo"): a direcao truncada nos m primeiros termos, como o FORM com M termos;
      o CoV por linha contra m e contra a fracao de alpha^2 que a direcao ve;
   6. dependencia do tempo: a falha e o menor g ao longo de 0..48 h; uma rodada de line sampling da a curva P_f(T)
      inteira. *)
SetDirectory[DirectoryName[$InputFileName]];
Quiet[CreateDirectory["figs_ls"]];
phi[x_] := CDF[NormalDistribution[], x];
cauda[c_] := phi[-c];                                   (* P(N(0,1) > c) = Phi(-c) *)

(* ---------- parte 1: a figura do metodo ---------- *)
th = 55 Degree; beta = 2.0; kap = 0.13;
al = {Cos[th], Sin[th]}; tg = {-Sin[th], Cos[th]};       (* e = -alpha = xi*/beta (para a falha) e o plano ortogonal *)
ponto[s_, r_] := s tg + r al;
cDe[s_] := beta + kap s^2;                              (* a superficie de falha: r = beta + kap s^2 *)
zs = MapThread[ponto, {{-2.3, -1.25, -0.35, 0.55, 1.45, 2.3}, {0.7, -0.9, 1.1, -0.35, 0.45, -1.3}}];
ss = zs . tg; cs = cDe /@ ss; i1 = 5;
falha = RegionPlot[{x, y} . al > cDe[{x, y} . tg], {x, -4, 4}, {y, -4, 5},
   PlotStyle -> Directive[Lighter[Red, 0.85]], BoundaryStyle -> None, PlotPoints -> 60];
superficie = ParametricPlot[ponto[s, cDe[s]], {s, -4, 4}, PlotStyle -> Directive[Red, Thick]];
esquerda = Show[falha, superficie,
   Graphics[{
     {GrayLevel[0.8], Circle[{0, 0}, #] & /@ {1, 2, 3}},
     {Black, Dashed, Line[{ponto[-4.5, 0], ponto[4.5, 0]}]},
     Table[With[{zp = ss[[i]] tg, xi = ponto[ss[[i]], cs[[i]]]},
       {{GrayLevel[0.55], Thin, Line[{zp - 1.2 al, xi}]},
        {Red, Thickness[0.006], Line[{xi, zp + 5 al}]},
        {Blue, Dotted, Line[{zs[[i]], zp}]},
        {Blue, PointSize[0.018], Point[zs[[i]]]},
        {White, EdgeForm[Blue], Disk[zp, 0.07]},
        {Red, PointSize[0.02], Point[xi]}}], {i, Length[zs]}],
     {Darker[Green], Thickness[0.006], Arrowheads[0.04], Arrow[{{0, 0}, beta al}]},
     {Black, PointSize[0.022], Point[{0, 0}], Point[beta al]},
     Text[Style["\!\(\*SuperscriptBox[\(\[Xi]\), \(*\)]\) = \[Beta] e", 13, Darker[Green], Background -> White], beta al - 0.75 tg - 0.35 al],
     Text[Style["e", 14, Darker[Green], Bold], 0.9 al + 0.25 tg],
     Text[Style["origem", 11], {0, 0} - 0.3 al + 0.35 tg],
     Text[Style["plano ortogonal a e", 11, Italic], ponto[-2.6, -0.3]],
     Text[Style["falha\ng(\[Xi]) \[LessEqual] 0", 13, Darker[Red], Background -> White], {2.75, 3.75}],
     Text[Style["\!\(\*SubscriptBox[\(\[Xi]\), \(\[ScriptL]\)]\) ~ N(0, I)", 12, Blue], zs[[i1]] + {0.2, 0.3}],
     Text[Style["\!\(\*SubsuperscriptBox[\(\[Xi]\), \(\[ScriptL]\), \(\[UpTee]\)]\)", 13, Blue], ss[[i1]] tg - 0.35 al - 0.1 tg],
     Text[Style["\!\(\*SubscriptBox[\(\[Nu]\), \(\[ScriptL]\)]\)", 14, Red, Background -> White], ponto[ss[[i1]], cs[[i1]]] + 0.4 tg]}],
   PlotRange -> {{-3.6, 3.6}, {-3.2, 4.2}}, AspectRatio -> Automatic, Frame -> True,
   FrameLabel -> {"\!\(\*SubscriptBox[\(\[Xi]\), \(1\)]\)", "\!\(\*SubscriptBox[\(\[Xi]\), \(2\)]\)"},
   PlotLabel -> Style["(a) cada amostra vira uma linha paralela a e", 13], ImageSize -> 470];
direita = Show[
   Plot[PDF[NormalDistribution[], c], {c, -3.5, 4.5}, PlotStyle -> Directive[Black, Thick], PlotRange -> {0, 0.45}],
   Plot[PDF[NormalDistribution[], c], {c, cs[[i1]], 4.5}, PlotStyle -> None, Filling -> Axis, FillingStyle -> Lighter[Red, 0.5]],
   Graphics[{
     {Red, Thick, Line[{{cs[[i1]], 0}, {cs[[i1]], 0.33}}]},
     {Darker[Green], Dashed, Line[{{beta, 0}, {beta, 0.4}}]},
     Text[Style["\!\(\*SubscriptBox[\(\[Nu]\), \(\[ScriptL]\)]\)", 14, Red], {cs[[i1]] + 0.2, 0.35}],
     Text[Style["\[Beta]  (FORM)", 13, Darker[Green]], {beta - 0.45, 0.42}],
     Text[Style["\!\(\*SubscriptBox[\(p\), \(\[ScriptL]\)]\) = \[CapitalPhi](-\!\(\*SubscriptBox[\(\[Nu]\), \(\[ScriptL]\)]\))", 14, Darker[Red]], {cs[[i1]] + 1.0, 0.1}],
     Text[Style["ao longo da linha, \[Nu] ~ N(0, 1)", 12, Italic], {-1.2, 0.43}],
     Text[Style["\!\(\*SubscriptBox[OverscriptBox[\(P\), \(^\)], \(f\)]\) = \!\(\*FractionBox[\(1\), \(N\)]\) \!\(\*UnderoverscriptBox[\(\[Sum]\), \(\[ScriptL] = 1\), \(N\)]\) \[CapitalPhi](-\!\(\*SubscriptBox[\(\[Nu]\), \(\[ScriptL]\)]\))", 16], {-1.6, 0.25}],
     Text[Style["s\[OAcute] a varia\[CCedilla]\[ATilde]o dos \!\(\*SubscriptBox[\(\[Nu]\), \(\[ScriptL]\)]\)\nentre as linhas entra na vari\[AHat]ncia", 11], {-1.6, 0.11}]}],
   Frame -> True, FrameLabel -> {"\[Nu] (coordenada ao longo de e)", "densidade"},
   PlotLabel -> Style["(b) cada linha d\[AAcute] a sua probabilidade exata", 13], ImageSize -> 430, AspectRatio -> 1.05];
Export["figs_ls/ls_esquema.pdf", GraphicsRow[{esquerda, direita}, Spacings -> 10, ImageSize -> 920]];
Print["parte 1: figs_ls/ls_esquema.pdf"];

(* ---------- parte 2: superficie curva, LS x AI x MC, mesmo numero de analises ---------- *)
(* falha: u_n > beta + kap y^2, com u_n = xi . e e y = xi . t (e = (0, 1), t = (1, 0)); beta = 2.5, kap = 0.15 *)
b2 = 2.5; k2 = 0.15;
pfExata2 = NIntegrate[PDF[NormalDistribution[], y] cauda[b2 + k2 y^2], {y, -Infinity, Infinity}];
Print["parte 2: P_f exata = ", N[pfExata2, 6], ";  FORM Phi(-beta) = ", N[cauda[b2], 6]];
nAn = 200; nRep = 200;                                 (* 200 analises por estimativa, 200 repeticoes *)
SeedRandom[1];
lsRep = Table[Mean[cauda[b2 + k2 #^2] & /@ RandomVariate[NormalDistribution[], nAn/2]], {nRep}];   (* 2 analises por linha *)
aiRep = Table[With[{x = RandomVariate[NormalDistribution[], {nAn, 2}] + ConstantArray[{0, b2}, nAn]},
    Mean[Boole[#[[2]] > b2 + k2 #[[1]]^2] Exp[-b2 #[[2]] + b2^2/2] & /@ x]], {nRep}];   (* AI N(xi*, I) *)
mcRep = Table[With[{x = RandomVariate[NormalDistribution[], {nAn, 2}]}, Mean[Boole[#[[2]] > b2 + k2 #[[1]]^2] & /@ x]], {nRep}];
tab2 = Table[{nome, Mean[r]/pfExata2 - 1, StandardDeviation[r]/Mean[r]}, {nome, {"line sampling", "AI no ponto de projeto", "MC direto"}},
   {r, {Switch[nome, "line sampling", lsRep, "AI no ponto de projeto", aiRep, "MC direto", mcRep]}}][[All, 1]];
Print["parte 2 (", nAn, " analises, ", nRep, " repeticoes): {metodo, erro da media, CoV} = ", N[tab2, 3]];
(* a convergencia de uma rodada de cada, para a figura *)
SeedRandom[7];
ys = RandomVariate[NormalDistribution[], 500];
pls = cauda[b2 + k2 #^2] & /@ ys;
xai = RandomVariate[NormalDistribution[], {1000, 2}] + ConstantArray[{0, b2}, 1000];
iai = Boole[#[[2]] > b2 + k2 #[[1]]^2] Exp[-b2 #[[2]] + b2^2/2] & /@ xai;
acum[v_, custo_] := Transpose[{custo Range[Length[v]], Accumulate[v]/Range[Length[v]]}];
covAcum[v_, custo_] := Table[{custo k, StandardDeviation[v[[;; k]]]/(Sqrt[k] Mean[v[[;; k]]])}, {k, 5, Length[v]}];
gConv = ListLogLinearPlot[{acum[pls, 2], acum[iai, 1]}, Joined -> True, PlotStyle -> {Blue, Orange},
   PlotRange -> {All, {0, 2.5 pfExata2}}, GridLines -> {None, {pfExata2}}, GridLinesStyle -> Directive[Black, Dashed],
   Frame -> True, FrameLabel -> {"an\[AAcute]lises (g)", "\!\(\*SubscriptBox[OverscriptBox[\(P\), \(^\)], \(f\)]\)"},
   PlotLegends -> Placed[{"line sampling (2 por linha)", "AI no ponto de projeto"}, {0.7, 0.85}],
   PlotLabel -> Style["(a) superf\[IAcute]cie curva: a estimativa (tracejado: exata)", 12], ImageSize -> 430];
gCov = ListLogLogPlot[{covAcum[pls, 2], covAcum[iai, 1]}, Joined -> True, PlotStyle -> {Blue, Orange},
   GridLines -> {None, {0.05}}, GridLinesStyle -> Directive[Black, Dashed], PlotRange -> {All, {0.01, 1}},
   Frame -> True, FrameLabel -> {"an\[AAcute]lises (g)", "CoV estimado"},
   PlotLegends -> Placed[{"line sampling", "AI no ponto de projeto"}, {0.7, 0.85}],
   PlotLabel -> Style["(b) o CoV (tracejado: 5%)", 12], ImageSize -> 430];

(* ---------- parte 3: superficie plana, variancia zero ---------- *)
(* g = b2 - xi2 (normal n = (0, 1)). Linhas sorteadas com a direcao e girada de d graus: na reta zp + nu e, a raiz e
   nu = (b2 - zp2)/e2. Com e = n (d = 0), todas as linhas dao nu = b2: variancia zero. Girada, a media continua exata,
   mas as linhas se espalham *)
lsPlanaRod[d_, nl_] := Module[{e = {Sin[d Degree], Cos[d Degree]}, z, zp},
   Table[z = RandomVariate[NormalDistribution[], 2]; zp = z - (z . e) e; cauda[(b2 - zp[[2]])/e[[2]]], {nl}]];
SeedRandom[2];
p0 = Flatten[Table[lsPlanaRod[0, 50], {5}]];
p20 = Table[lsPlanaRod[20, 50], {200}];
Print["parte 3: plana, e = normal, 5 rodadas de 50 linhas: menor e maior p_l = ", N[{Min[p0], Max[p0]}, 8], " (exata ",
  N[cauda[b2], 8], "), desvio entre as linhas = ", StandardDeviation[p0], ";  e girada de 20 graus, 200 rodadas: erro da media = ",
  N[Mean[Mean /@ p20]/cauda[b2] - 1, 3], ", CoV por linha = ", N[StandardDeviation[Flatten[p20]]/Mean[Flatten[p20]], 3]];

(* ---------- parte 4: direcao errada, sem vies ---------- *)
(* a direcao do line sampling girada de d graus em relacao a normal. Na reta zp + nu e, a falha
   (zp2 + nu e2) > b2 + k2 (zp1 + nu e1)^2 e uma quadratica concava em nu: a linha esta na falha entre as duas raizes
   nu1 < nu2 (com e1 = 0, uma raiz so, e p = Phi(-nu1)). Entao p = Phi(-nu1) - Phi(-nu2), e 0 sem raiz real: com a direcao
   girada, a linha pode entrar e sair da falha *)
lsGirada[d_, nl_] := Module[{e = {Sin[d Degree], Cos[d Degree]}, z, zp, qa, qb, qc, disc, r},
   Table[z = RandomVariate[NormalDistribution[], 2]; zp = z - (z . e) e;
    qa = -k2 e[[1]]^2; qb = e[[2]] - 2 k2 zp[[1]] e[[1]]; qc = zp[[2]] - b2 - k2 zp[[1]]^2;
    If[qa == 0, cauda[-qc/qb],
     disc = qb^2 - 4 qa qc;
     If[disc <= 0, 0., r = Sort[{(-qb - Sqrt[disc])/(2 qa), (-qb + Sqrt[disc])/(2 qa)}]; cauda[r[[1]]] - cauda[r[[2]]]]], {nl}]];
SeedRandom[3];
tab4 = Table[With[{r = Table[Mean[lsGirada[d, 50]], {100}]}, {d, Mean[r]/pfExata2 - 1, StandardDeviation[r]/Mean[r] Sqrt[50]}],
   {d, {0, 10, 20, 30}}];
Print["parte 4: {graus, erro da media de 100 rodadas de 50 linhas, CoV por linha} = ", N[tab4, 3]];

(* ---------- parte 5: muitas dimensoes, a direcao truncada ---------- *)
(* g = b0 - a . xi, com n = 2000 e a_j ~ Sqrt[lambda_j] de um espectro que decai (como os modos de um campo).
   P_f exata = Phi(-b0/|a|). A direcao do FORM com m termos e e_m = a_(1..m)/|a_(1..m)|: ao longo dela, g cai com
   inclinacao |a_(1..m)|, e o resto de a fica no plano ortogonal e faz cada linha ter o seu nu. *)
n5 = 2000; lam5 = Table[Exp[-j/60.] + 0.02 Exp[-j/600.], {j, n5}]; a5 = Sqrt[lam5]; a5 = 1.0 a5/Norm[a5];   (* |a| = 1 *)
b5 = 2.5; pf5 = cauda[b5];
covLinha5[m_] := Module[{am = Join[a5[[;; m]], ConstantArray[0., n5 - m]], e, sa, srest},
   e = am/Norm[am]; sa = a5 . e;                         (* inclinacao ao longo de e *)
   srest = Sqrt[1 - sa^2];                               (* |parte de a fora de e| *)
   (* nu_l = (b0 - a . zp)/sa, com a . zp ~ N(0, srest^2): p_l = Phi(-nu_l); CoV por linha exato por integracao *)
   With[{m1 = NIntegrate[PDF[NormalDistribution[], v] cauda[(b5 - srest v)/sa], {v, -Infinity, Infinity}],
      m2 = NIntegrate[PDF[NormalDistribution[], v] cauda[(b5 - srest v)/sa]^2, {v, -Infinity, Infinity}]},
    {m, sa^2, m1/pf5 - 1, Sqrt[m2 - m1^2]/m1}]];
tab5 = covLinha5 /@ {5, 10, 20, 50, 100, 200, 500, 1000, 2000};
Print["parte 5: {m, alfa^2 visto pela direcao, erro da media (exato), CoV por linha (exato)} = ", N[tab5, 3]];
gTrunc = ListLogLogPlot[{Most[tab5][[All, {1, 4}]]}, Joined -> True, PlotMarkers -> Automatic, PlotStyle -> Blue,
   Frame -> True, FrameLabel -> {"m (termos da dire\[CCedilla]\[ATilde]o do FORM)", "CoV por linha"},
   Epilog -> {Text[Style[ToString[Round[100 #[[2]]]] <> "%", 9, Gray], {Log[#[[1]]], Log[#[[4]]]}, {-0.3, -1.6}] & /@ Most[tab5]},
   PlotLabel -> Style["(c) n = 2000, dire\[CCedilla]\[ATilde]o truncada em m termos\n(r\[OAcute]tulos: a fra\[CCedilla]\[ATilde]o de |a\!\(\*SuperscriptBox[\(|\), \(2\)]\) que ela v\[EHat]; com m = 2000, CoV = 0)", 11],
   ImageSize -> 430];

(* ---------- parte 6: dependencia do tempo ---------- *)
(* g(xi, t) = beta(t) - (cos th(t) xi1 + sin th(t) xi2), t = 0..48 h de hora em hora: dois episodios criticos (12 h e
   36 h), com mecanismos em direcoes diferentes. Falha ate T: min_{t <= T} g < 0. Line sampling com a direcao do
   instante mais critico (36 h); cada linha da nu_l(T) para todo T com as mesmas analises (o g(t) de todas as horas).
   Aqui cada hora e linear em nu, entao a raiz de cada hora sai das mesmas 2 analises da linha. *)
tt = Range[0, 48]; bT[t_] := 3.0 - 0.55 Exp[-((t - 12)/5.)^2] - 0.7 Exp[-((t - 36)/5.)^2];
thT[t_] := (15 + 70/(1 + Exp[-(t - 24)/2.])) Degree;     (* dois episodios: mecanismo a 15 graus ate ~24 h, a 85 depois *)
gTempo[x_] := Table[bT[t] - (Cos[thT[t]] x[[1]] + Sin[thT[t]] x[[2]]), {t, tt}];
e6 = {Cos[thT[36]], Sin[thT[36]]};                      (* a direcao do instante mais critico, 36 h *) t6 = {-e6[[2]], e6[[1]]};
SeedRandom[5];
nl6 = 1000;
res6 = Table[Module[{z = RandomVariate[NormalDistribution[], 2], zp, gc},
    zp = z - (z . e6) e6;
    (* ao longo da linha, cada hora e linear em nu: g_t(nu) = g_t(zp) - (e . d_t) nu; a raiz da hora t *)
    gc = Table[With[{d = {Cos[thT[t]], Sin[thT[t]]}}, (bT[t] - d . zp)/(d . e6)], {t, tt}];
    {FoldList[Min, gc], gc}], {nl6}];                    (* nu_l(T) = a menor raiz ate T, e a raiz de cada hora *)
cT = res6[[All, 1]];
pfInstLS = Table[Mean[cauda /@ res6[[All, 2, k]]], {k, Length[tt]}];   (* a instantanea pelas mesmas linhas *)
pfT = Table[Mean[cauda /@ cT[[All, k]]], {k, Length[tt]}];
cvT = Table[StandardDeviation[cauda /@ cT[[All, k]]]/(Sqrt[nl6] pfT[[k]]), {k, Length[tt]}];
SeedRandom[6];
xmc = RandomVariate[NormalDistribution[], {200000, 2}];
gmc = gTempo /@ xmc;
pfMC = Table[Mean[N[Boole[Min[#[[;; k]]] < 0] & /@ gmc]], {k, Length[tt]}];
pfInst = Table[cauda[bT[t]], {t, tt}];
Print["parte 6: instantanea maxima pelas linhas = ", N[Max[pfInstLS], 4], ";  P_f(48 h) acumulada: LS (", nl6, " linhas, ", 2 nl6, " analises) = ", N[pfT[[-1]], 4], " (CoV ", N[cvT[[-1]], 2],
  ");  MC de 200000 = ", N[pfMC[[-1]], 4], ";  P_f instantanea maxima (36 h) = ", N[Max[pfInst], 4]];
gTempoFig = ListLogPlot[{Transpose[{tt, pfMC}], Transpose[{tt, pfT}], Transpose[{tt, pfInst}], Transpose[{tt, pfInstLS}][[;; ;; 2]]},
   Joined -> {True, False, True, False},
   PlotStyle -> {Directive[Black, Thick], Blue, Directive[Gray, Dashed], Darker[Red]}, PlotMarkers -> {None, {Automatic, 6}, None, {"\[EmptyCircle]", 12}},
   PlotRange -> {All, {1*^-4, 0.1}}, Frame -> True, FrameLabel -> {"T (h)", "probabilidade de falha"},
   PlotLegends -> Placed[{"acumulada at\[EAcute] T, MC de 200 000", "acumulada, line sampling (1000 linhas)", "instant\[AHat]nea, exata", "instant\[AHat]nea, line sampling"}, {0.65, 0.22}],
   PlotLabel -> Style["(d) com o tempo: a curva \!\(\*SubscriptBox[\(P\), \(f\)]\)(T) inteira sai das mesmas linhas", 12], ImageSize -> 430];
Export["figs_ls/ls_exemplos.pdf", GraphicsGrid[{{gConv, gCov}, {gTrunc, gTempoFig}}, Spacings -> {20, 20}, ImageSize -> 900]];
Print["figuras: figs_ls/ls_esquema.pdf e figs_ls/ls_exemplos.pdf"];
