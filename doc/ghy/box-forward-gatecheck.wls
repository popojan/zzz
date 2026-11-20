(* Gate-check: does box-averaging (per-term sinc weights) improve method B's
   zero LOCATION at kappa = gap*log X in (pi, 3pi), and is it null at kappa <= pi?
   Hypothesis (stated before running):
     H1: at kappa ~ 2pi-3pi, median |t*-gamma| improves for c ~ 0.25-0.5 (w = c*gap).
     H2: at kappa ~ pi, gain ~ 1x (band saturation: in-band reweighting is inert).
   Box-average of F_B over [t-w,t+w] == multiply each p^m term by Sinc[w*m*Log[p]].
   Target stays k-1/2: box-avg of the true staircase crosses k-1/2 exactly at
   gamma_k when the window holds one zero (bias only from close pairs).
   N0 box-average correction ~ w^2/(12 pi t) ~ 1e-6 here: ignored. *)

kLo = 1000; kHi = 1100;
zs = Table[N[Im[ZetaZero[k]]], {k, kLo - 1, kHi + 1}];
gm[k_] := zs[[k - (kLo - 1) + 1]];
gapAt[t_] := 2 Pi/Log[t/(2 Pi)];
n0[t_] := t/(2 Pi) Log[t/(2 Pi E)] + 7/8;

ppData[X_] := Module[{ps = Prime[Range[PrimePi[X]]]},
  Transpose[Flatten[
    Table[{m Log[p], N[p^(-m/2)/m]}, {p, ps}, {m, 1, Floor[Log[X]/Log[p]]}], 1]]];

locate[k_, freqs_, amps_, c_] := Module[
  {g = gm[k], lo, hi, w, wt, fb, flo, fhi, mid, fm},
  lo = 0.5 (gm[k - 1] + gm[k]); hi = 0.5 (gm[k] + gm[k + 1]);
  w = c gapAt[g];
  wt = If[c == 0, ConstantArray[1., Length[freqs]], Sinc[w freqs]];
  fb[t_] := n0[t] - (wt amps).Sin[t freqs]/Pi - (k - 0.5);
  flo = fb[lo]; fhi = fb[hi];
  If[flo fhi > 0,
    lo = lo - 0.5 gapAt[g]; hi = hi + 0.5 gapAt[g];
    flo = fb[lo]; fhi = fb[hi]];
  If[flo fhi > 0, Return[$Failed]];
  Do[mid = 0.5 (lo + hi); fm = fb[mid];
     If[flo fm <= 0, hi = mid; fhi = fm, lo = mid; flo = fm], {60}];
  Abs[0.5 (lo + hi) - g]/gapAt[g]];

cs = {0, 0.125, 0.25, 0.375, 0.5, 0.75};
Do[
  Module[{fr, am, kap, res, meds, fails},
    {fr, am} = ppData[X];
    kap = Log[X] gapAt[gm[(kLo + kHi)/2]];
    Print["X=", X, "  #terms=", Length[fr], "  kappa=", NumberForm[kap/Pi, {4, 2}], " pi"];
    Do[
      res = Table[locate[k, fr, am, c], {k, kLo, kHi}];
      fails = Count[res, $Failed];
      meds = Median[DeleteCases[res, $Failed]];
      Print["  c=", NumberForm[c, {3, 3}], "  median|err|/gap=",
        NumberForm[meds, {6, 4}], "  RMS=",
        NumberForm[Sqrt[Mean[DeleteCases[res, $Failed]^2]], {6, 4}],
        "  fails=", fails],
      {c, cs}]],
  {X, {15, 59, 229, 3000}}]
