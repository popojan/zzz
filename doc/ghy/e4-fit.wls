(* e4-fit.wls -- Step 3 of the Weil-formula program: the in-band fit.
   (doc/notes/superres-multi-cutoff.md)

   Mechanism: method B's location-error field is band-limited (it is the
   truncation deficit, frequencies <= a = log X), and the zeros sample it
   below Nyquist exactly when gap*a < pi. Hence B's error is IN-BAND and
   recoverable by direct inversion on the well-conditioned channels of
   exact Weil functionals -- no fixed point, no soft modes, no kernel
   model error. The evanescent (sub-Rayleigh) component stays at a
   Tikhonov prior anchored to the B seeds.

   Functionals, single centre Tc (exact rational), even/odd completions:
     h_cos(r): bumps W(x)cos(wx) at r = +-Tc  ->  zero side 2 Sum W(x_j)cos(w x_j)
     h_sin(r): odd-completed                   ->  zero side 2 Sum W(x_j)sin(w x_j)
   with x_j = gamma_j - Tc, W(x) = sinc(b x)^(2qw), band delta = 2 qw b.
   Weil pairs (gHat = B-spline of order 2qw, support delta):
     g_cos(u) = cos(Tc u) [gHat(u-w) + gHat(u+w)]
     g_sin(u) = sin(Tc u) [gHat(u+w) - gHat(u-w)]
   admissible iff w + delta <= a. Archimedean term at height Tc ~ 1e21
   (psi-asymptotics, x-Taylor; corrections O(1/Tc) ~ 1e-21 dropped):
     arch_cos = 2 log(Tc/(2 pi)) gHat(w),   arch_sin = 0.
   Pole terms ~ W(Tc) ~ 0. So the DATA identity reads
     2 Sum_j W(x_j) cos/sin(w x_j) = arch - 2 Sum_{p^m<=X} (log p) p^(-m/2) g(u),
   computable from primes alone; only cos/sin(Tc u) needs high precision
   (phase folding); everything else runs in machine arithmetic.

   Part 0 validates the g_cos/g_sin algebra as an identity at T = 143.5
   against zeros2 (1000-digit zeros) with the numeric-psi archimedean term.
   Part 1 fits window zeros at 10^22+5000 (k = 10^4) from B seeds
   (/tmp/bseeds_e22k1e4.tsv) and scores against Odlyzko zeros5.
*)

prec = 100;
kP = 10000; X = Prime[kP]; aBand = Log[X];
qw = 8; delta = 3; bw = delta/(2 qw);

parseDec[s_String] := Module[{t, parts, ip, fp},
  t = StringTrim[s]; parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1, FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  ip + fp];
numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];

(* B-spline pair *)
bspline[m_, y_] := If[Abs[y] >= m/2, 0,
  Sum[(-1)^k Binomial[m, k] With[{t = y + m/2 - k}, If[t > 0, t^(m - 1), 0]],
    {k, 0, m}]/(m - 1)!];
WW[x_] := If[x == 0, 1, (Sin[bw x]/(bw x))^(2 qw)];
dWW[x_] := If[x == 0, 0,
  2 qw (Sin[bw x]/(bw x))^(2 qw - 1) (bw Cos[bw x]/(bw x) - bw Sin[bw x]/(bw x)^2)];
gHat[u_] := bspline[2 qw, u/(2 bw)]/(2 bw);

(* prime-power tables: u = m log p (high precision), amp = (log p) p^(-m/2) *)
pps = Flatten[Table[
    Table[{m N[Log[p], prec + 30], N[Log[p] p^(-m/2), 40]}, {m, 1, Floor[Log[X]/Log[p]]}],
    {p, Prime[Range[kP]]}], 1];
uHi = pps[[All, 1]]; amp = pps[[All, 2]];
Print["prime powers: ", Length[pps]];

(* data D_alpha for centre Tc (exact rational) and frequency grid ws:
   phases folded at high precision once, the rest machine *)
weilData[Tc_, ws_] := Module[{ph, cT, sT, um, gp, gm, dcos, dsin, arch},
  ph = Mod[N[Tc, prec + 30] uHi, 2 Pi];
  cT = N[Cos[ph]]; sT = N[Sin[ph]];
  um = N[uHi];
  arch = 2 Log[N[Tc]/(2 Pi)];
  Table[
   gp = gHat /@ (um - w); gm = gHat /@ (um + w);
   dcos = arch gHat[w] - 2 (amp cT) . (gp + gm);
   dsin = 0 - 2 (amp sT) . (gm - gp);
   {dcos, dsin}, {w, ws}]];

(* ============ Part 0: identity check at low height ============ *)
(* light settings: the check only needs ~10 digits to catch algebra bugs;
   skip with `wolframscript -file e4-fit.wls part1` *)
runPart0 = !MemberQ[Rest[$ScriptCommandLine], "part1"];
If[runPart0,
Print["=== Part 0: family identity check, T=143.5, X=10^5"];
Module[{lines, zstrs = {}, cur = "", g2, Tc = 1435/10, ws = {0, 3, 7}, dat,
    zside, psiInt, diffs},
  lines = ReadList[FileNameJoin[{scriptDir, "..", "refs", "odlyzko", "zeros2"}], String, NullRecords -> True];
  Do[Module[{l = StringTrim[line]},
    Which[l === "", If[cur =!= "", AppendTo[zstrs, cur]; cur = ""],
      StringContainsQ[l, "."], If[cur =!= "", AppendTo[zstrs, cur]]; cur = l,
      True, cur = cur <> l]], {line, lines}];
  If[cur =!= "", AppendTo[zstrs, cur]];
  g2 = parseDec /@ zstrs;
  dat = weilData[Tc, ws];
  (* low height: numeric psi instead of the asymptotic arch; subtract the
     asymptotic arch already inside weilData and add the exact integral *)
  diffs = Table[Module[{w = ws[[i]], zc, zsn, archAsym, archCos, archSin, xm},
     zc = 2 Sum[WW[N[g - Tc, 60]] Cos[w N[g - Tc, 60]], {g, g2}];
     zsn = 2 Sum[WW[N[g - Tc, 60]] Sin[w N[g - Tc, 60]], {g, g2}];
     archAsym = 2 Log[N[Tc]/(2 Pi)] gHat[w];
     xm = 50;
     archCos = (1/Pi) NIntegrate[WW[x] Cos[w x] (Re[PolyGamma[0, 1/4 + I (Tc + x)/2]] - Log[Pi]),
        Evaluate[Flatten[{x, Range[-xm, xm, 5]}]],
        WorkingPrecision -> 30, PrecisionGoal -> 12, AccuracyGoal -> 14, MaxRecursion -> 8];
     archSin = (1/Pi) NIntegrate[WW[x] Sin[w x] (Re[PolyGamma[0, 1/4 + I (Tc + x)/2]] - Log[Pi]),
        Evaluate[Flatten[{x, Range[-xm, xm, 5]}]],
        WorkingPrecision -> 30, PrecisionGoal -> 12, AccuracyGoal -> 14, MaxRecursion -> 8];
     {zc - (dat[[i, 1]] - archAsym + archCos), zsn - (dat[[i, 2]] + archSin)}],
    {i, Length[ws]}];
  Print["  identity residuals (cos, sin) per w ", ws, ":"];
  Print["  ", ScientificForm[N[diffs], 3]]]];

(* ============ Part 1: the fit at 10^22 + 5000 ============ *)
Print["=== Part 1: in-band fit at 10^22+5000, k=10^4, gap*a = 1.55"];
b22 = 1370919909931995300000;
z5 = parseDec /@ Select[ReadList[FileNameJoin[{scriptDir, "..", "refs", "odlyzko", "zeros5"}], String], numLineQ];
bs = Association[(ToExpression[#[[1]]] -> parseDec[#[[2]]]) & /@
   (StringSplit[#, "\t"] & /@ ReadList["/tmp/bseeds_e22k1e4.tsv", String])];
offs = Sort[Keys[bs]];
Print["B seeds: offsets ", Min[offs], "..", Max[offs]];

center = 5000;
(* bs values are ABSOLUTE ordinates (zzz output); z5 are offsets from b22 *)
Tc = Round[bs[center] 10^6]/10^6;   (* exact rational centre *)
xB = Association[Table[o -> N[bs[o] - Tc, 30], {o, offs}]];          (* B offsets *)
xT = Association[Table[o -> N[b22 + z5[[o]] - Tc, 30], {o, offs}]];  (* truth     *)

winW = 12; farW = 20;   (* B-error leakage through W: ~4e-4 at +-8, ~4e-7 at +-12 *)
fitO = Select[offs, Abs[xB[#]] <= winW &];
farO = Select[offs, winW < Abs[xB[#]] <= farW &];
M = Length[fitO];
Print["fit zeros M=", M, "  far fixed: ", Length[farO]];

wmax = aBand - delta - 1/20;
ws = Range[0, N[wmax], N[wmax]/56];
dat = weilData[Tc, ws];
Print[Length[ws], " frequencies, ", 2 Length[ws], " functionals"];

(* model and Jacobian, machine arithmetic *)
farCos = Table[2 Sum[WW[xB[o]] Cos[w xB[o]], {o, farO}], {w, ws}];
farSin = Table[2 Sum[WW[xB[o]] Sin[w xB[o]], {o, farO}], {w, ws}];
model[xs_] := Flatten[Table[
   {farCos[[i]] + 2 Sum[WW[xs[[j]]] Cos[ws[[i]] xs[[j]]], {j, M}],
    farSin[[i]] + 2 Sum[WW[xs[[j]]] Sin[ws[[i]] xs[[j]]], {j, M}]}, {i, Length[ws]}]];
jac[xs_] := Module[{rows},
  rows = Table[Flatten[Table[
      {2 (dWW[xs[[j]]] Cos[ws[[i]] xs[[j]]] - WW[xs[[j]]] ws[[i]] Sin[ws[[i]] xs[[j]]]),
       2 (dWW[xs[[j]]] Sin[ws[[i]] xs[[j]]] + WW[xs[[j]]] ws[[i]] Cos[ws[[i]] xs[[j]]])},
      {i, Length[ws]}]], {j, M}];
  Transpose[rows]];

dataVec = Flatten[dat];
sigD = 10.^-6;        (* far-zone + arithmetic floor *)
sigP = 0.012;         (* B prior *)
Off[LinearSolve::luc];
x0 = N[xB /@ fitO];
xs = x0;
(* normal equations rescaled by sigD^2 and solved at 50 digits: the ridge
   (sigD/sigP)^2 ~ 7e-9 against rank-deficient J^T J is beyond machine *)
Do[Module[{r, J, A, rhs, dx},
   r = dataVec - model[xs];
   J = jac[xs];
   A = Transpose[J] . J + (sigD/sigP)^2 IdentityMatrix[M];
   rhs = Transpose[J] . r - (sigD/sigP)^2 (xs - x0);
   dx = N[LinearSolve[SetPrecision[A, 50], SetPrecision[rhs, 50]]];
   xs = xs + dx;
   Print["  GN iter ", it, ": |r| = ", ScientificForm[Norm[dataVec - model[xs]], 3],
     "  max|dx| = ", ScientificForm[Max[Abs[dx]], 3]]],
  {it, 16}];

(* score: central zone |x| <= 4 *)
ctr = Select[Range[M], Abs[x0[[#]]] <= 4 &];
errB = Table[xB[fitO[[j]]] - xT[fitO[[j]]], {j, ctr}];
errF = Table[xs[[j]] - xT[fitO[[j]]], {j, ctr}];
Print["central zone (", Length[ctr], " zeros):"];
Print["  B    mean|err| = ", ScientificForm[Mean[Abs[errB]], 4],
  "   max = ", ScientificForm[Max[Abs[errB]], 3]];
Print["  fit  mean|err| = ", ScientificForm[Mean[Abs[errF]], 4],
  "   max = ", ScientificForm[Max[Abs[errF]], 3]];
Print["  gain = ", N[Mean[Abs[errB]]/Mean[Abs[errF]], 4],
  "x   (true-seed ceiling 3.5x, B abs 0.0119)"];

(* control: does the data distinguish truth from the fit/B at all? *)
xsT = N[xT /@ fitO];
Print["control |r|:  at fit = ", ScientificForm[Norm[dataVec - model[xs]], 4],
  "   at TRUTH = ", ScientificForm[Norm[dataVec - model[xsT]], 4],
  "   at B = ", ScientificForm[Norm[dataVec - model[x0]], 4]];
Print["corr(fit err, B err) central = ", N[Correlation[errF, errB], 4]];
