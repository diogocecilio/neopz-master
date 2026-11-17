(* ::Package:: *)
(* exemplo_1d.wl: a amostragem por importancia (AI) e a variavel de controle num exemplo pequeno, a partir do codigo do
   Diogo: X ~ N(6, 1), falha quando X < 1. Cada parte remete a uma equacao da formulacao (formulacao.pdf, Apendice A).
   Rodar:  wolframscript -file exemplo_1d.wl > exemplo_1d_saida.txt   (ou colar parte por parte num notebook) *)

SeedRandom[2026];
fi[x_] := PDF[NormalDistribution[0, 1], x];   (* phi, a densidade normal padrao *)
Fi[x_] := CDF[NormalDistribution[0, 1], x];   (* Phi, a acumulada *)
sci[x_, d_: 3] := If[x == 0, "0", Module[{e = Floor[Log10[Abs[x]]]},
    ToString[NumberForm[N[x/10^e], {d + 2, d}]] <> "e" <> ToString[e]]];   (* 2.770e-7 *)
fix[x_, d_: 1] := ToString[NumberForm[N[x], {12, d}]];                       (* 5906.4 *)

(* ======== Parte 1: o codigo original, a AI no espaco X ======== *)
pdf[x_] := PDF[NormalDistribution[6, 1], x];   (* f: a densidade original de X *)
pdf2[x_] := PDF[NormalDistribution[2, 1], x];  (* h: a densidade de amostragem *)
n = 10000;                                      (* inteiro: RandomVariate pede 10000, e nao 10000. *)
hsamples = RandomVariate[NormalDistribution[2, 1], n];
fail = 0;
w = {};
For[i = 1, i <= Length[hsamples], i++,
  If[hsamples[[i]] < 1,
    hx = pdf2[hsamples[[i]]];
    fx = pdf[hsamples[[i]]];
    AppendTo[w, fx/hx];
    fail++;
  ];
];
pAI = Total[w]/Length[hsamples];
pExata = NIntegrate[pdf[x], {x, -Infinity, 1}];
Print["Parte 1. AI no espaco X: ", fail, " falhas em ", n, " amostras; P_f = ", sci[pAI], "; exata = ", sci[pExata],
  " = Phi(-5) = ", sci[Fi[-5]]];

(* ======== Parte 2: o mesmo no espaco padrao (o Y do Beck, o xi da formulacao) ========
   xi = X - 6 ~ N(0, 1). A falha X < 1 vira xi < -5: beta = 5, e o ponto de projeto e xi* = -5.
   Sortear X de N(2, 1) e o mesmo que sortear xi de N(mu, 1), com mu = 2 - 6 = -4. *)
mu = -4.;
(* 2a. A linha do "xi = mu + eps": mu + eps, com eps ~ N(0, 1), tem distribuicao N(mu, 1) *)
Print["Parte 2a. a lei de mm + eps, eps ~ N(0, 1): ",
  TransformedDistribution[mm + e, e \[Distributed] NormalDistribution[0, 1]]];
(* 2b. A Eq. (28) em 1D: o peso e a razao entre a densidade original e a de amostragem, no mesmo ponto s *)
Print["Parte 2b. f(s)/h(s) = ", Simplify[PDF[NormalDistribution[0, 1], s]/PDF[NormalDistribution[mm, 1], s]]];
(* 2c. As duas contas dao o mesmo peso nas mesmas amostras *)
xis = hsamples - 6;                         (* as mesmas amostras, no espaco padrao *)
wX = pdf[hsamples]/pdf2[hsamples];          (* f/h no espaco X, como na Parte 1 *)
wXi = Exp[-mu xis + mu^2/2];                (* a Eq. (28) *)
Print["Parte 2c. maior diferenca relativa entre f/h e exp(-mu xi + mu^2/2): ", sci[Max[Abs[wX/wXi - 1]], 1]];
(* 2d. Com xi = mu + eps: h(xi) = phi(eps), mas f(xi) = phi(mu + eps), que nao e phi(eps). Por isso o peso nao e 1 *)
eps = xis - mu;
Print["Parte 2d. h(xi) = phi(eps)? ", Max[Abs[PDF[NormalDistribution[mu, 1], xis] - fi[eps]]] < 10^-15,
  ";  f(xi) = phi(eps)? ", Max[Abs[fi[xis] - fi[eps]]] < 10^-15,
  ";  w = exp(-mu eps - mu^2/2)? ", Max[Abs[wXi/Exp[-mu eps - mu^2/2] - 1]] < 10^-12];
Print["           P_f pela Eq. (27) no espaco padrao: ", sci[Mean[UnitStep[-5 - xis] wXi]], " (a mesma da Parte 1)"];

(* ======== Parte 3: qualquer mu da a resposta certa (sem vies); o mu so muda a variancia ========
   A variancia exata de uma amostra (Eq. 27): E_h[I w^2] - P_f^2, com E_h[I w^2] = Exp[mu^2] Phi(mu - beta). *)
beta = 5.;
pf = Fi[-beta];
ai[m_, nn_] := Module[{x = RandomVariate[NormalDistribution[m, 1], nn]}, Mean[UnitStep[-beta - x] Exp[-m x + m^2/2]]];
covExato[m_, nn_] := Sqrt[(Exp[m^2] Fi[m - beta] - pf^2)/nn]/pf;
Print["Parte 3. 200 repeticoes de 1000 amostras, P_f exata = ", sci[pf]];
Do[est = Table[ai[m, 1000], {200}];
  Print["  mu = ", fix[m], ": media ", sci[Mean[est]], ",  CoV observado ", fix[100 StandardDeviation[est]/pf],
    "%,  exato ", fix[100 covExato[m, 1000]], "%"],
  {m, {0., -2., -3., -4., -5., -6., -7.}}];
muCE = -fi[beta]/Fi[-beta];                    (* E[xi | falha]: o mu da entropia cruzada (Secao 5) *)
muMin = t /. Last[FindMinimum[covExato[t, 1000], {t, -5.}]];
Print["  E[xi | falha] = ", fix[muCE, 3], ";  o mu de menor variancia = ", fix[muMin, 3], ";  o ponto de projeto = -5"];

(* ======== Parte 4: a variavel de controle ========
   Para ela ter o que fazer, o "FS" e uma funcao nao linear de duas normais padrao: xi, que a AI desloca e que o
   substituto ve, e zeta, que faz o papel dos modos que o substituto nao ve (fica N(0, 1), sem peso, como os modos que
   o FORM nao desloca na Secao 5.2):
       ln FS(xi, zeta) = 0,04 (xi + 5) + 0,004 (xi + 5)^2 + kz zeta,   kz = 0,002.
   FS ~ 1,35 na media; a falha, FS < 1, fica perto de xi = -5. O kz mede o que o substituto nao ve: com kz = 0,002 o
   ganho e ~2 (rho ~ 0,77, como no talude); com kz = 0,01 ele cai para ~1, porque perto da fronteira, onde os pesos sao
   maiores, quem decide a falha passa a ser o zeta. *)
kz = 0.002;
lnFS[x_, z_] := 0.04 (x + 5) + 0.004 (x + 5)^2 + kz z;
pf4 = NIntegrate[Fi[-(0.04 (x + 5) + 0.004 (x + 5)^2)/kz] fi[x], {x, -10, 10}];
Print["Parte 4. P_f exata (integral em xi, com zeta integrado a mao) = ", sci[pf4]];
mu4 = -5.;                                      (* o deslocamento: o ponto de projeto de xi *)
umaAI[nn_] := Module[{xp, zp, lfp, b0, b1, mP, tau, xs, zs, w4, iF, ghat, h, EH, y, c, gam},
  (* 4a. o substituto linear de ln FS (Eq. 34), por minimos quadrados num piloto separado de 200 amostras da AI:
     assim ele nao usa as amostras do estimador (e o papel da validacao cruzada no codigo) *)
  xp = RandomVariate[NormalDistribution[mu4, 1], 200]; zp = RandomVariate[NormalDistribution[0, 1], 200];
  lfp = MapThread[lnFS, {xp, zp}];
  {b0, b1} = LinearModelFit[Transpose[{xp, lfp}], s, s]["BestFitParameters"];
  ghat[x_] := b0 + b1 x;
  (* 4b. o indicador suavizado (Eq. 35): ln FS = ghat + residuo, residuo ~ N(m, tau^2), entao
     P(falha | xi) = Phi((-ghat - m)/tau). Aqui m e tau sao a media e o desvio dos residuos do ajuste; no codigo
     (diagnostico_ai.py) saem por maxima verossimilhanca nas falhas (o probit) *)
  {mP, tau} = {Mean[#], StandardDeviation[#]} &[lfp - ghat[xp]];
  h[x_] := Fi[(-ghat[x] - mP)/tau];
  EH = Fi[(-b0 - mP)/Sqrt[tau^2 + b1^2]];       (* a media exata de H sob phi (Eq. 36) *)
  (* 4c. a AI: so xi e deslocado; o peso so tem xi (Eq. 28) *)
  xs = RandomVariate[NormalDistribution[mu4, 1], nn]; zs = RandomVariate[NormalDistribution[0, 1], nn];
  w4 = Exp[-mu4 xs + mu4^2/2];
  iF = UnitStep[-MapThread[lnFS, {xs, zs}]];
  (* 4d. o estimador com controle (Eq. 33): I w e H w, com gamma da regressao de I w em H w *)
  y = iF w4; c = h[xs] w4;
  gam = Covariance[y, c]/Variance[c];
  {Mean[y], Mean[y - gam c] + gam EH, StandardDeviation[y]/Sqrt[nn], StandardDeviation[y - gam c]/Sqrt[nn],
   Correlation[y, c], EH, NIntegrate[h[x] fi[x], {x, -Infinity, Infinity}]}];
r1 = umaAI[1000];
Print["  uma AI de 1000: P = ", sci[r1[[1]]], " +- ", sci[r1[[3]], 2], ";  com controle ", sci[r1[[2]]], " +- ",
  sci[r1[[4]], 2], ";  rho = ", fix[r1[[5]], 2]];
Print["  E[H] pela Eq. (36) = ", sci[r1[[6]], 6], ";  por integracao numerica = ", sci[r1[[7]], 6]];
reps = Table[umaAI[1000][[1 ;; 5]], {200}];
Print["  200 AI de 1000 (cada uma com o seu piloto):"];
Print["    sem controle: media ", sci[Mean[reps[[All, 1]]]], ", desvio entre as AI ", sci[StandardDeviation[reps[[All, 1]]], 2],
  " (estimado ", sci[Mean[reps[[All, 3]]], 2], ")"];
Print["    com controle: media ", sci[Mean[reps[[All, 2]]]], ", desvio entre as AI ", sci[StandardDeviation[reps[[All, 2]]], 2],
  " (estimado ", sci[Mean[reps[[All, 4]]], 2], ")"];
Print["    ganho = razao das variancias = ", fix[Variance[reps[[All, 1]]]/Variance[reps[[All, 2]]], 2],
  ";  1/(1 - rho^2) com o rho medio (Eq. 32) = ", fix[1/(1 - Mean[reps[[All, 5]]]^2), 2]];
