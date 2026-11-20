(* e3-multiL-fit.wls -- Experiment E3 of doc/notes/superres-multi-cutoff.md:
   multi-cutoff (multi-L) full-profile least squares at the 10^36 proxy regime
   (10^22, k <= 600). Extends e2-singleL-fit.wls: the residual runs over a 2D
   grid (T_i, L'_l); each cutoff gets its own band-edge nuisance block
   {cos(L'_l tau), sin(L'_l tau), tau cos(L'_l tau), tau sin(L'_l tau)}
   (far-zero leakage sits at the band edge, whose position moves with L'),
   eliminated by variable projection. A zero-position shift has a fixed,
   known Si((tau - g_j) L') signature across cutoffs -- the leverage that a
   single cutoff lacks (E2 result: gain 1.01x, i.e. none).

     model_l(tau) = Sum_j Si((tau - g_j) L'_l)/pi + nuisance_l(tau)
     J_l[i,j]     = -(L'_l/pi) sinc((tau_i - g_j) L'_l)

   Usage: wolframscript -file e3-multiL-fit.wls [Wfit] [halfSpan]
*)

prec = 50;
argv = Rest[$ScriptCommandLine];
Wfit = If[Length[argv] >= 1, ToExpression[argv[[1]]], 24];
halfSpan = If[Length[argv] >= 2, Rationalize[ToExpression[argv[[2]]], 0], 5/2];
ks = {150, 200, 270, 360, 450, 600};
b22 = 1370919909931995300000;
npts = 175;
sigmaT = 8/10;
targets = Range[50, 69];

Xs = Prime /@ ks; Lsm = N[Log[Xs]];
Print["Wfit=", Wfit, " halfSpan=", N[halfSpan], " cutoffs X=", Xs];

parseDec[s_String] := Module[{t, parts, ip, fp},
  t = StringTrim[s]; parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1, FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  ip + fp];
numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];
z5 = parseDec /@ Select[ReadList[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "..", "refs", "odlyzko", "zeros5"}], String], numLineQ];
truthAbs[o_] := b22 + z5[[o]];

bsRows = StringSplit[#, "\t"] & /@ ReadList["/tmp/bseeds_e22k600.tsv", String];
bseed = Association[(ToExpression[#[[1]]] -> parseDec[#[[2]]]) & /@ bsRows];

(* prime-power terms up to the largest cutoff; per-cutoff masked amplitudes *)
Xmax = Max[Xs];
terms = Flatten[Table[
    Table[{N[p^m], m N[Log[p], prec + 20], N[p^(-m/2)/m, prec]}, {m, 1, Floor[Log[Xmax]/Log[p]]}],
    {p, Prime[Range[Max[ks]]]}], 1];
mlogp = terms[[All, 2]];
ampsL = Table[If[terms[[t, 1]] <= Xs[[l]], terms[[t, 3]], 0], {l, Length[ks]}, {t, Length[terms]}];
Print["prime-power terms: ", Length[terms]];

n0smooth[T_] := With[{x = T/(2 Pi)}, x Log[x/E] + 7/8];

(* reduced data: rows stacked per cutoff, all O(1) in machine precision *)
reducedData[Tc_, taus_, n0m1_, M_] := Module[{Ts, phases, sins, im},
  Ts = Tc + taus;
  phases = Outer[Times, mlogp, N[Ts, prec + 20]];
  sins = Sin[phases];
  im = -(ampsL . sins);                       (* Length[ks] x npts *)
  N[Flatten[Table[
    n0smooth[N[Ts, prec + 10]] + im[[l]]/Pi - n0m1 - M/2, {l, Length[ks]}]]]];

si = SinIntegral;
Off[LinearSolve::luc];

fitOne[o_] := Module[{Tc, taus, n0, M, f, g, w, basis, r, J, lam, chi2, chi2new,
    gnew, r2, c2, dgs, nL = Length[ks], taum},
  M = 2 Wfit + 1;
  n0 = 10^22 + o - Wfit;
  Tc = bseed[o];
  taus = Table[-halfSpan + (i - 1) (2 halfSpan)/(npts - 1), {i, npts}];
  f = reducedData[Tc, taus, n0 - 1, M];
  taum = N[taus];
  g = N[Table[bseed[o + j] - Tc, {j, -Wfit, Wfit}]];
  w = Flatten[Table[Exp[-taum^2/(2 N[sigmaT]^2)], {nL}]];
  (* block-diagonal band-edge nuisance: 4 columns per cutoff *)
  basis = ArrayFlatten[Table[If[l == lc,
      Transpose[{Cos[Lsm[[l]] taum], Sin[Lsm[[l]] taum],
                 taum Cos[Lsm[[l]] taum], taum Sin[Lsm[[l]] taum]}],
      ConstantArray[0., {npts, 4}]], {l, nL}, {lc, nL}]];
  chi2Of[gv_] := Module[{m0, rr, cc},
    m0 = Flatten[Table[Sum[si[Lsm[[l]] (taum - gv[[j]])], {j, M}]/Pi, {l, nL}]];
    rr = f - m0;
    cc = LeastSquares[Sqrt[w] basis, Sqrt[w] rr];
    rr = rr - basis . cc;
    {w . (rr^2), rr, cc}];
  {chi2, r, c2} = chi2Of[g];
  lam = 10.^-3;
  Do[
    J = ArrayFlatten[Table[
      Table[-(Lsm[[l]]/Pi) Sinc[Lsm[[l]] (taum - g[[j]])], {j, M}] // Transpose, {l, nL}, {lc, 1}]];
    J = J - basis . LeastSquares[Sqrt[w] basis, Sqrt[w] J];
    dgs = LinearSolve[Transpose[J] . (w J) + lam IdentityMatrix[M], Transpose[J] . (w r)];
    gnew = g + dgs;
    {chi2new, r2, c2} = chi2Of[gnew];
    If[chi2new < chi2,
      g = gnew; {chi2, r} = {chi2new, r2}; lam = Max[lam/3, 10.^-8],
      lam = lam 10];
    If[Max[Abs[dgs]] < 10.^-7, Break[]],
    {25}];
  g[[Wfit + 1]]];

errsB = Table[N[bseed[o] - truthAbs[o]], {o, targets}];
Print["B baseline   mean|err| = ", Mean[Abs[errsB]], "  max = ", Max[Abs[errsB]]];
errsFit = Table[Module[{gf = fitOne[o], e},
   e = gf - N[truthAbs[o] - bseed[o]];
   Print["  o=", o, "  fit err = ", e, "   (B: ", N[bseed[o] - truthAbs[o]], ")"]; e], {o, targets}];
Print["E3 fit       mean|err| = ", Mean[Abs[errsFit]], "  max = ", Max[Abs[errsFit]]];
Print["gain vs B = ", Mean[Abs[errsB]]/Mean[Abs[errsFit]], "x   (true-seed bound: 2.7x)"];
