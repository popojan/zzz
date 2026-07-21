(* e5-above-threshold.wls -- the practical payoff experiment: exact
   Weil-identity fitting ABOVE the saturation threshold, where surplus
   information exists (band-saturation.md, "above-threshold constants").

   Regime: 10^12 + 5000, k = 10^5 (X = p_100000 = 1299709, a = log X ~ 14.08,
   kappa = gap*a ~ 3.69 > pi, ~17% DOF surplus). Here --boot already beats
   B, but it pays three losses the Weil fit does not: E_1-kernel model error
   (54-80% rms capture), crossing-reading (one number per zero), and
   truncated neighbour subtraction. Machinery identical to e4-fit.wls
   Part 1 (validated there + identity-checked at low height); only the
   configuration differs. Scored against Odlyzko zeros3 (floor 1e-8).

   Inputs: /tmp/bseeds_e12k1e5.tsv  (B locations, offsets 4850..5150, ABSOLUTE)
           /tmp/boot_e12k1e5.tsv    (boot 32 baseline, offsets 4995..5004)
           doc/refs/odlyzko/zeros3  (truth: line j = zero 10^12 + j, base below)
*)

prec = 100;
kP = If[Length[Rest[$ScriptCommandLine]] >= 1,
  ToExpression[Rest[$ScriptCommandLine][[1]]], 100000];
X = Prime[kP]; aBand = Log[X];
qw = 8; delta = 3; bw = delta/(2 qw);
b12 = 267653395647;

parseDec[s_String] := Module[{t, parts, ip, fp},
  t = StringTrim[s]; parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1, FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  ip + fp];
numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];

bspline[m_, y_] := If[Abs[y] >= m/2, 0,
  Sum[(-1)^k Binomial[m, k] With[{t = y + m/2 - k}, If[t > 0, t^(m - 1), 0]],
    {k, 0, m}]/(m - 1)!];
WW[x_] := If[x == 0, 1, (Sin[bw x]/(bw x))^(2 qw)];
dWW[x_] := If[x == 0, 0,
  2 qw (Sin[bw x]/(bw x))^(2 qw - 1) (bw Cos[bw x]/(bw x) - bw Sin[bw x]/(bw x)^2)];
gHat[u_] := bspline[2 qw, u/(2 bw)]/(2 bw);

Print["building prime tables, k=", kP, " X=", X, " a=", N[aBand, 6]];
pps = Flatten[Table[
    Table[{m N[Log[p], prec + 30], N[Log[p] p^(-m/2), 40]}, {m, 1, Floor[Log[X]/Log[p]]}],
    {p, Prime[Range[kP]]}], 1];
uHi = pps[[All, 1]]; amp = pps[[All, 2]];
Print["prime powers: ", Length[pps]];

weilData[Tc_, ws_] := Module[{ph, cT, sT, um, gp, gm, arch},
  ph = Mod[N[Tc, prec + 30] uHi, 2 Pi];
  cT = N[Cos[ph]]; sT = N[Sin[ph]];
  um = N[uHi];
  arch = 2 Log[N[Tc]/(2 Pi)];
  Table[
   gp = gHat /@ (um - w); gm = gHat /@ (um + w);
   {arch gHat[w] - 2 (amp cT) . (gp + gm),
    -2 (amp sT) . (gm - gp)}, {w, ws}]];

(* ---- configuration ---- *)
z3 = parseDec /@ Select[ReadList[FileNameJoin[{scriptDir, "..", "refs", "odlyzko", "zeros3"}], String], numLineQ];
bs = Association[(ToExpression[#[[1]]] -> parseDec[#[[2]]]) & /@
   (StringSplit[#, "\t"] & /@ ReadList["/tmp/bseeds_e12k1e5.tsv", String])];
bt = Association[(ToExpression[#[[1]]] -> parseDec[#[[2]]]) & /@
   (StringSplit[#, "\t"] & /@ ReadList["/tmp/boot_e12k1e5.tsv", String])];
offs = Sort[Keys[bs]];
Print["B seeds: ", Min[offs], "..", Max[offs], "   boot baseline: ", Length[bt], " zeros"];

center = 5000;
Tc = Round[bs[center] 10^6]/10^6;
xB = Association[Table[o -> N[bs[o] - Tc, 30], {o, offs}]];
xT = Association[Table[o -> N[b12 + z3[[o]] - Tc, 30], {o, offs}]];

winW = 12; farW = 20;
fitO = Select[offs, Abs[xB[#]] <= winW &];
farO = Select[offs, winW < Abs[xB[#]] <= farW &];
M = Length[fitO];
Print["fit zeros M=", M, "  far fixed: ", Length[farO], "  kappa = ",
  N[2 Pi/Log[N[Tc]/(2 Pi)] aBand, 4]];

wmax = aBand - delta - 1/20;
ws = Range[0, N[wmax], N[wmax]/56];
dat = weilData[Tc, ws];
Print[Length[ws], " frequencies, ", 2 Length[ws], " functionals"];

farCos = Table[2 Sum[WW[xB[o]] Cos[w xB[o]], {o, farO}], {w, ws}];
farSin = Table[2 Sum[WW[xB[o]] Sin[w xB[o]], {o, farO}], {w, ws}];
model[xs_] := Flatten[Table[
   {farCos[[i]] + 2 Sum[WW[xs[[j]]] Cos[ws[[i]] xs[[j]]], {j, M}],
    farSin[[i]] + 2 Sum[WW[xs[[j]]] Sin[ws[[i]] xs[[j]]], {j, M}]}, {i, Length[ws]}]];
jac[xs_] := Transpose[Table[Flatten[Table[
     {2 (dWW[xs[[j]]] Cos[ws[[i]] xs[[j]]] - WW[xs[[j]]] ws[[i]] Sin[ws[[i]] xs[[j]]]),
      2 (dWW[xs[[j]]] Sin[ws[[i]] xs[[j]]] + WW[xs[[j]]] ws[[i]] Cos[ws[[i]] xs[[j]]])},
     {i, Length[ws]}]], {j, M}]];

dataVec = Flatten[dat];
sigD = 10.^-7;
sigP = 0.006;
Off[LinearSolve::luc];
x0 = N[xB /@ fitO];
xs = x0;
(* Levenberg-Marquardt with step clipping: higher cutoffs sharpen the
   functionals (linearization radius ~ 1/w_max), undamped GN diverges *)
lam = 10.^-3; stepMax = 0.3 N[2 Pi/Log[N[Tc]/(2 Pi)]];
Do[Module[{r, J, JtJ, A, rhs, dx, xsNew, rNew, tries = 0},
   r = dataVec - model[xs];
   J = jac[xs];
   JtJ = Transpose[J] . J;
   While[tries < 8,
    A = JtJ + (sigD/sigP)^2 IdentityMatrix[M] + lam DiagonalMatrix[Diagonal[JtJ] + 10.^-12];
    rhs = Transpose[J] . r - (sigD/sigP)^2 (xs - x0);
    dx = N[LinearSolve[SetPrecision[A, 50], SetPrecision[rhs, 50]]];
    If[Max[Abs[dx]] > stepMax, dx = dx stepMax/Max[Abs[dx]]];
    xsNew = xs + dx;
    rNew = dataVec - model[xsNew];
    If[Norm[rNew] < Norm[r], xs = xsNew; lam = Max[lam/3, 10.^-7]; Break[],
      lam *= 10; ++tries]];
   Print["  LM iter ", it, ": |r| = ", ScientificForm[Norm[dataVec - model[xs]], 3],
     "  max|dx| = ", ScientificForm[Max[Abs[dx]], 3], "  lam = ", ScientificForm[lam, 2]]],
  {it, 24}];

(* ---- scoring ---- *)
score[zone_] := Module[{idx, eB, eF},
  idx = Select[Range[M], Abs[x0[[#]]] <= zone &];
  eB = Table[xB[fitO[[j]]] - xT[fitO[[j]]], {j, idx}];
  eF = Table[xs[[j]] - xT[fitO[[j]]], {j, idx}];
  Print["zone |x|<=", zone, "  (", Length[idx], " zeros):"];
  Print["  B    mean|err| = ", ScientificForm[N[Mean[Abs[eB]]], 4],
    "   max = ", ScientificForm[N[Max[Abs[eB]]], 3]];
  Print["  fit  mean|err| = ", ScientificForm[N[Mean[Abs[eF]]], 4],
    "   max = ", ScientificForm[N[Max[Abs[eF]]], 3],
    "   gain vs B = ", NumberForm[N[Mean[Abs[eB]]/Mean[Abs[eF]]], 4], "x"]];
score[6]; score[3];

(* three-way on the boot-baseline zeros *)
bo = Sort[Keys[bt]];
Print["three-way on offsets ", Min[bo], "..", Max[bo], ":"];
eB3 = Table[N[bs[o] - (b12 + z3[[o]])], {o, bo}];
eC3 = Table[N[bt[o] - (b12 + z3[[o]])], {o, bo}];
eF3 = Table[Module[{j = Position[fitO, o][[1, 1]]}, N[xs[[j]] - xT[o]]], {o, bo}];
Print["  B    mean|err| = ", ScientificForm[Mean[Abs[eB3]], 4]];
Print["  boot mean|err| = ", ScientificForm[Mean[Abs[eC3]], 4],
  "   (", NumberForm[Mean[Abs[eB3]]/Mean[Abs[eC3]], 4], "x over B)"];
Print["  fit  mean|err| = ", ScientificForm[Mean[Abs[eF3]], 4],
  "   (", NumberForm[Mean[Abs[eB3]]/Mean[Abs[eF3]], 4], "x over B, ",
  NumberForm[Mean[Abs[eC3]]/Mean[Abs[eF3]], 4], "x over boot)"];
Print["  per-zero |err| (B / boot / fit):"];
Do[Module[{j = Position[fitO, bo[[i]]][[1, 1]]},
  Print["   o=", bo[[i]], "  ", ScientificForm[Abs[eB3[[i]]], 3], "  ",
    ScientificForm[Abs[eC3[[i]]], 3], "  ", ScientificForm[Abs[eF3[[i]]], 3]]],
  {i, Length[bo]}];
