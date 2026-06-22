#!/usr/bin/env wolframscript
(* CAN A BETTER BACKWARD DETECTOR RAISE alpha ON NOISY ZEROS?

   reach-exponent.wls measured alpha ~ 0.74 on EXACT zeros (super-critical) but
   the loop on its own noisy self-zeros sits ~0.5 (marginal).  The gap is the
   tax of zero-approximation noise + a naive peak-threshold detector.  This
   script tests whether sidelobe-suppressing windows (the snap-margin lever,
   which extended REACH under noise at fixed n) raise the reach EXPONENT alpha
   on the loop's ACTUAL noisy zeros (/tmp/loop60k.txt -- the real method-B error
   field, not synthetic).

   Detector: psi'_w(x)/env, window w(gamma) on the cos-sum, DC-normalized
   (Total[w]=N) so peak height matches across windows; threshold tau=0.30.
   Windows: Dirichlet (=naive), Fejer (triangle, non-neg kernel), Jackson
   (Fejer^2, stronger sidelobe kill), Gauss.  Reach = flawless prefix (no missed
   prime, no composite false-positive; higher prime powers don't-care), scored
   against PrimeQ/PrimePowerQ (scoring only).  alpha = slope of log reach vs log N.  *)

dir = DirectoryName[$InputFileName];
gAll = N[ToExpression /@ StringTrim /@ ReadList["/tmp/loop60k.txt", String]];

env[x_, n_] := (20/4) (LogIntegral[n - 4] Log[x])/x;
tau = 0.30; Xcap = 2200; blk = 300;
norm[w_, n_] := w n/Total[w];

windowVec[name_, gv_] := Module[{g1 = gv/Last[gv], n = Length[gv]},
  norm[Switch[name,
     "Dirichlet", ConstantArray[1., n],
     "Fejer",     1 - g1,
     "Jackson",   (1 - g1)^2,
     "Gauss",     Exp[-(2.0 g1)^2]], n]];

(* flawless-prefix reach for the first n loop zeros with window `name` *)
reach[n_, name_] := Module[{gv = Take[gAll, n], w, xm = 1},
  w = windowVec[name, gv];
  Catch[Do[
    Module[{xs = Range[x0, Min[x0 + blk - 1, Xcap]], lx, M, peak, r, isP, isPP},
     lx = Log[N[xs]];
     M = Cos[Outer[Times, lx, gv]];
     peak = 1 + 1/(xs - xs^3) - (4/Sqrt[N[xs]]) (M . w);
     r = peak/env[N[xs], n];
     isP = PrimeQ /@ xs; isPP = PrimePowerQ /@ xs;
     Do[If[(isP[[i]] && r[[i]] <= tau) || (! isPP[[i]] && r[[i]] > tau),
           Throw[xs[[i]] - 1]]; xm = xs[[i]], {i, Length[xs]}]],
    {x0, 2, Xcap, blk}]; xm]];

ngrid = {750, 1500, 3000, 6000};
names = {"Dirichlet", "Fejer", "Jackson"};
data = AssociationMap[Function[nm, reach[#, nm] & /@ ngrid], names];
WriteString["stdout", "# (progress) reaches computed\n"];

WriteString["stdout", "# BACKWARD REACH on the loop's OWN noisy zeros (gamma<=" <>
   ToString[NumberForm[Last[gAll], 6]] <> ", tau=0.30).  reach = flawless prefix.\n\n"];
WriteString["stdout", "n       " <> StringJoin[StringPadLeft[#, 11] & /@ names] <> "\n"];
Do[WriteString["stdout", StringPadRight[ToString[ngrid[[i]]], 8] <>
   StringJoin[StringPadLeft[ToString[data[#][[i]]] <>
       If[data[#][[i]] >= Xcap, "(cap)", ""], 11] & /@ names] <> "\n"], {i, Length[ngrid]}];

WriteString["stdout", "\n# reach exponent alpha (log reach vs log n, global fit)\n"];
Do[Module[{a = LinearModelFit[
     Transpose[{Log[N@ngrid], Log[N@data[name]]}], xx, xx]["BestFitParameters"][[2]]},
   WriteString["stdout", StringPadRight[name, 11] <> "alpha = " <>
      ToString[NumberForm[a, {4, 3}]] <> "\n"]], {name, names}];

WriteString["stdout", "\n# alpha > 0.5 => loop super-critical at kappa=pi on its OWN noisy zeros.\n"];
WriteString["stdout", "# Best window's alpha vs Dirichlet's = the recoverable backward lever.\n"];

s = OpenWrite[FileNameJoin[{dir, "reach-detectors.tsv"}]];
WriteString[s, "n\t" <> StringRiffle[names, "\t"] <> "\n"];
Do[WriteString[s, ToString[ngrid[[i]]] <> "\t" <>
   StringRiffle[ToString[data[#][[i]]] & /@ names, "\t"] <> "\n"], {i, Length[ngrid]}];
Close[s];
