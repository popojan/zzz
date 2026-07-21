#!/usr/bin/env wolframscript
(* BACKWARD REACH EXPONENT -- is the zeros->primes step information-limited
   (a real wall, like the forward step) or merely DETECTOR-limited (slack)?

   The loop's growth rate is X_{i+1} ~ X_i^{alpha * 2pi/kappa}, where alpha is
   the backward reach exponent: from n zeros (height T), the simple psi'
   threshold detector reaches primes up to x_max ~ n^alpha.  At kappa=pi the
   loop is super-critical iff alpha > 1/2.  Memory says alpha ~ 0.5-0.6 (naive
   detector).  But the EXPLICIT FORMULA determines primes up to x ~ T (Nyquist:
   truncation error of psi(x) = x - sum_{|gamma|<T} x^rho/rho is ~ x (log x)^2/T,
   < 1 for x < T/(log T)^2) -- nearly LINEAR in n, far above n^0.5.  So if the
   measured alpha is ~0.5 while the information ceiling is ~1, the backward step
   is DETECTOR-limited and a better detector could push the loop super-critical.

   This measures alpha for the project's actual (envelope-threshold) detector on
   EXACT zeros (no approximation noise -- isolates the detector vs the
   information ceiling), and reports reach / gamma_n (Nyquist fraction).  *)

dir   = DirectoryName[$InputFileName];
zfile = FileNameJoin[{dir, "..", "refs", "odlyzko", "zeros1"}];
allg  = N[ToExpression /@ StringTrim /@ ReadList[zfile, String]];

env[x_, n_] := (20/4) (LogIntegral[n - 4] Log[x])/x;
tau = 0.30;
Xcap = 6000;                        (* scan ceiling *)
blk  = 400;                         (* x-block size (vectorized) *)

(* flawless-prefix reach for n exact zeros: scan x upward IN BLOCKS, stop just
   before the first MISSED PRIME (prime, peak<=tau) or FALSE POSITIVE (genuine
   composite, peak>tau).  higher prime powers are don't-cares (loop bookkeeping). *)
reach[n_] := Module[{g = Take[allg, n], xm = 1},
  Catch[Do[
    Module[{xs = Range[x0, Min[x0 + blk - 1, Xcap]], lx, M, peak, r, isP, isPP},
     lx = Log[N[xs]];
     M = Cos[Outer[Times, lx, g]];                       (* |xs| x n *)
     peak = 1 + 1/(xs - xs^3) - (4/Sqrt[N[xs]]) (M . ConstantArray[1., n]);
     r = peak/env[N[xs], n];
     isP = PrimeQ /@ xs; isPP = PrimePowerQ /@ xs;
     Do[If[(isP[[i]] && r[[i]] <= tau) || (! isPP[[i]] && r[[i]] > tau),
           Throw[xs[[i]] - 1]]; xm = xs[[i]], {i, Length[xs]}]],
    {x0, 2, Xcap, blk}]; xm]];

ngrid = {500, 1000, 2000, 4000, 8000, 16000, 32000};
res = {#, allg[[#]], reach[#]} & /@ ngrid;

WriteString["stdout", "# BACKWARD REACH EXPONENT  (exact zeros, envelope detector, tau=0.30)\n"];
WriteString["stdout", "#   information ceiling (explicit-formula Nyquist): reach ~ gamma_n  (frac -> O(1))\n\n"];
WriteString["stdout", "n        gamma_n     reach x_max   x_max/gamma_n   alpha(local)\n"];
Do[Module[{n = res[[i, 1]], gn = res[[i, 2]], xm = res[[i, 3]], a = ""},
   If[i > 1, a = ToString[NumberForm[
        N[Log[xm/res[[i - 1, 3]]]/Log[n/res[[i - 1, 1]]]], {4, 3}]]];
   WriteString["stdout",
     StringPadRight[ToString[n], 9] <>
     StringPadRight[ToString[NumberForm[gn, {7, 1}]], 12] <>
     StringPadRight[ToString[xm] <> If[xm >= Xcap, "(cap)", ""], 14] <>
     StringPadRight[ToString[NumberForm[N[xm/gn], {5, 4}]], 16] <>
     a <> "\n"]], {i, Length[res]}];

(* global fit alpha over the whole range *)
fit = LinearModelFit[Table[{Log[res[[i, 1]]], Log[res[[i, 3]]]}, {i, Length[res]}], xx, xx];
WriteString["stdout", "\n# global fit:  x_max ~ n^alpha,  alpha = " <>
   ToString[NumberForm[fit["BestFitParameters"][[2]], {4, 3}]] <> "\n"];
WriteString["stdout", "# alpha > 0.5  => loop super-critical at kappa=pi (frontier shifts unboundedly).\n"];
WriteString["stdout", "# x_max/gamma_n << 1  => detector-limited, NOT information-limited: slack remains.\n"];
