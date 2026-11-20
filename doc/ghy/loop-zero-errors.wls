#!/usr/bin/env wolframscript
(* Distribution of the gap-normalized error of the self-computed zeros in a
   `zzz --loop` checkpoint, against exact ZetaZero (one call per zero -- hence
   sampled).  normErr_k = (approx_k - Im rho_k) / gap_k, gap_k = 2 pi / log(gamma_k/2 pi)
   -- i.e. the location error in units of the local mean spacing.

     wolframscript -file doc/ghy/loop-zero-errors.wls [state] [n_samples]
   Writes loop-zero-errors.png (histogram + error-vs-height scatter).            *)

file = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "zzz-loop.state"];
K    = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 600];

lines = ReadList[file, String];
zeros = ToExpression /@ Select[lines,                       (* zero lines: digits + a '.' *)
   StringMatchQ[#, (DigitCharacter | ".") ..] && StringContainsQ[#, "."] &];
n = Length[zeros];
If[n < 10, Print["too few zeros in ", file]; Exit[1]];

idx = DeleteDuplicates@Round@Subdivide[1, n, Min[K, n] - 1];   (* uniform sample of indices *)
Print["sampling ", Length[idx], " of ", n, " zeros (exact ZetaZero) ..."];

exact = Table[Im[ZetaZero[k]], {k, idx}];
approx = zeros[[idx]];
gap = 2 Pi/Log[exact/(2 Pi)];
nerr = (approx - exact)/gap;                                   (* normalized error *)

Print["normalized error  (units of local mean gap):"];
Print["  mean   = ", N@Mean[nerr]];
Print["  median = ", N@Median[nerr]];
Print["  stddev = ", N@StandardDeviation[nerr]];
Print["  max|.| = ", N@Max[Abs[nerr]]];
Print["absolute error |approx-exact|:  median = ", N@Median[Abs[approx - exact]],
      "   max = ", N@Max[Abs[approx - exact]]];

out = "loop-zero-errors.png";
img = GraphicsRow[{
   Histogram[nerr, Automatic, "PDF",
     Frame -> True, FrameLabel -> {"(approx - exact)/gap", "density"},
     PlotLabel -> Row[{"normalized zero error  (", Length[idx], " of ", n, ")"}]],
   ListPlot[Transpose[{exact, nerr}], Frame -> True,
     FrameLabel -> {"height \[Gamma]", "(approx - exact)/gap"},
     PlotLabel -> "error vs height", PlotRange -> All, PlotStyle -> PointSize[0.005]]
   }, ImageSize -> 1100];
Export[out, img];
Print["wrote ", out];
