(* e2-singleL-fit.wls -- Experiment E2 of doc/notes/superres-multi-cutoff.md:
   single-L full-profile least squares at the 10^36 proxy regime (10^22, k=600).

   Model (window-relative coordinates tau = T - Tc, g_j = gamma_j - Tc):
     f(tau)   = F(T, L) - (n0 - 1) - M/2
     model(g) = Sum_j [ Si((tau - g_j) L) / pi + 1/2 sign-shift folded into M/2 ]
                + alpha cos(L tau) + beta sin(L tau)        (band-edge nuisance)
   J_ij = d model_i / d g_j = -(L/pi) sinc((tau_i - g_j) L)
   Data from exact prime sums at 50 digits; LM on machine-precision reduced field.

   Inputs: /tmp/bseeds_e22k600.tsv (B locations, offsets 26..93, from zzz --ghy)
           doc/refs/odlyzko/zeros5 (truth, scoring only)
   Output: per-target fitted error vs B baseline; summary stats.
*)

prec = 50;
kPrimes = 600; X = Prime[kPrimes]; L = N[Log[X], prec]; Lm = N[L];
b22 = 1370919909931995300000;
Wfit = 24;            (* fitted zeros each side of target: M = 2 Wfit + 1 *)
halfSpan = 5/2;       (* T-grid half-width *)
npts = 175;
sigmaT = 8/10;        (* Gaussian weight width *)
targets = Range[50, 69];

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

(* prime-power terms: {m log p, p^(-m/2)/m} for p^m <= X *)
terms = Flatten[Table[
    Table[{m N[Log[p], prec + 20], N[p^(-m/2)/m, prec]}, {m, 1, Floor[Log[X]/Log[p]]}],
    {p, Prime[Range[kPrimes]]}], 1];
mlogp = terms[[All, 1]]; amp = terms[[All, 2]];
Print["prime-power terms: ", Length[terms]];

n0smooth[T_] := With[{x = T/(2 Pi)}, x Log[x/E] + 7/8];

(* reduced data field on the grid, exact-precision then machine *)
reducedData[Tc_, taus_, n0m1_, M_] := Module[{Ts, phases, imlogp, Fs},
  Ts = Tc + taus;
  phases = Outer[Times, mlogp, N[Ts, prec + 20]];
  imlogp = -(amp . Sin[phases]);
  Fs = n0smooth[N[Ts, prec + 10]] + imlogp/Pi;
  N[Fs - n0m1 - M/2]];

si = SinIntegral;

fitOne[o_] := Module[{Tc, taus, n0, M, f, g, w, basis, r, J, lam, chi2, chi2new,
    gnew, iter, A, c, r2, c2, dgs, Tcrat},
  M = 2 Wfit + 1;
  n0 = 10^22 + o - Wfit;            (* ordinal of first window zero *)
  Tcrat = Rationalize[bseed[o], 0]; (* grid center = B seed of target: no truth used *)
  Tc = Tcrat;
  taus = Table[-halfSpan + (i - 1) (2 halfSpan)/(npts - 1), {i, npts}];
  f = reducedData[Tc, taus, n0 - 1, M];
  taus = N[taus];
  g = N[Table[bseed[o + j] - Tc, {j, -Wfit, Wfit}]];
  w = Exp[-taus^2/(2 N[sigmaT]^2)];
  basis = Transpose[{Cos[Lm taus], Sin[Lm taus]}];
  chi2Of[gv_] := Module[{m0, rr, cc},
    m0 = Sum[si[Lm (taus - gv[[j]])], {j, Length[gv]}]/Pi;
    rr = f - m0;
    cc = LeastSquares[Sqrt[w] basis, Sqrt[w] rr];
    rr = rr - basis . cc;
    {w . (rr^2), rr, cc}];
  {chi2, r, c} = chi2Of[g];
  lam = 10.^-3;
  Do[
    J = Table[-(Lm/Pi) Sinc[Lm (taus - g[[j]])], {j, M}] // Transpose;
    (* project nuisance out of Jacobian as well *)
    J = J - basis . LeastSquares[Sqrt[w] basis, Sqrt[w] J];
    A = Transpose[J] . (w J);
    dgs = LinearSolve[A + lam IdentityMatrix[M], Transpose[J] . (w r)];
    gnew = g + dgs;
    {chi2new, r2, c2} = chi2Of[gnew];
    If[chi2new < chi2,
      g = gnew; {chi2, r, c} = {chi2new, r2, c2}; lam = Max[lam/3, 10.^-8],
      lam = lam 10];
    If[Max[Abs[dgs]] < 10.^-7, Break[]],
    {iter, 25}];
  g[[Wfit + 1]]];   (* window-relative: adding Tc back would drown O(0.01) in machine eps at 10^21 *)

Off[LinearSolve::luc];   (* expected: super-resolution conditioning; LM damping handles it *)

errsB = Table[N[bseed[o] - truthAbs[o]], {o, targets}];
Print["B baseline   mean|err| = ", Mean[Abs[errsB]], "  max = ", Max[Abs[errsB]]];
errsFit = Table[Module[{gf = fitOne[o], e},
   e = gf - N[truthAbs[o] - bseed[o]];   (* both terms O(1), exact difference first *)
   Print["  o=", o, "  fit err = ", e, "   (B: ", N[bseed[o] - truthAbs[o]], ")"]; e], {o, targets}];
Print["E2 fit       mean|err| = ", Mean[Abs[errsFit]], "  max = ", Max[Abs[errsFit]]];
Print["gain vs B = ", Mean[Abs[errsB]]/Mean[Abs[errsFit]], "x   (true-seed bound: 2.7x)"];
