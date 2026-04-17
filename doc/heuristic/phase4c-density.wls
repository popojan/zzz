#!/usr/bin/env wolframscript

(* ================================================================ *)
(* Phase 4c: Error analysis with density normalization              *)
(* Key question: is there a k-error relationship when we account   *)
(* for the zero density at height T?                               *)
(* ================================================================ *)

fmt[x_, n_:5] := ToString[AccountingForm[N[x], n]];

Print[""];
Print["Phase 4c: Density-Normalized Error Analysis"];
Print["============================================="];
Print[""];

nSmooth[t_] := 1 + RiemannSiegelTheta[t]/Pi;

primeSum[t_, k_] := Module[{s = 0.0 + 0.0 I},
  Do[s += Log[1 - Prime[i]^(-1/2 + I t)], {i, k}];
  -Im[s]/Pi
];

(* Zero spacing at height T *)
spacing[t_] := 2 Pi / Log[t/(2 Pi)];

(* ================================================================ *)
(* 1: Error in S(T) units and in spacing units                     *)
(* ================================================================ *)

Print["--- 1: Error vs k at various heights ---"];
Print[""];
Print["Columns: absolute S(T) error, and error as fraction of zero spacing"];
Print["(error/spacing < 0.5 means we locate zeros within half a spacing)"];
Print[""];

heights = {100, 500, 1000, 5000, 10000, 50000, 100000};
kVals = {10, 50, 100, 200, 500, 1000};

(* Header *)
Print[StringJoin[
  StringPadRight["", 22],
  StringJoin[StringPadRight["k=" <> ToString[#], 14] & /@ kVals]
]];
Print[StringJoin[Table["-", 106]]];

allData = {};

Do[
  t1 = N[Im[ZetaZero[nz]]];
  t2 = N[Im[ZetaZero[nz + 1]]];
  tmid = (t1 + t2)/2;
  sExact = nz - nSmooth[tmid];
  sp = spacing[tmid];

  row1 = "n=" <> StringPadRight[ToString[nz], 7] <>
    "T~" <> StringPadRight[fmt[tmid, 5], 10] <> "  ";
  row2 = StringPadRight["", 22];

  vals1 = {};
  vals2 = {};
  Do[
    fk = primeSum[tmid, k];
    err = Abs[fk - sExact];
    errNorm = err / sp;
    AppendTo[vals1, StringPadRight[fmt[err, 3], 14]];
    AppendTo[vals2, StringPadRight[fmt[errNorm, 3], 14]];
    AppendTo[allData, {nz, tmid, k, err, errNorm, sp}];
  , {k, kVals}];

  Print[row1, StringJoin[vals1], "  |err|"];
  Print[row2, StringJoin[vals2], "  err/spacing"];
  Print[""];
, {nz, heights}];

(* ================================================================ *)
(* 2: Focus on the error/spacing metric                            *)
(* ================================================================ *)

Print["--- 2: Mean error/spacing vs k ---"];
Print[""];

Do[
  subset = Select[allData, #[[3]] == k &];
  norms = subset[[All, 5]];
  Print["  k=", StringPadRight[ToString[k], 5],
    "  mean(err/spacing)=", fmt[Mean[norms], 4],
    "  max=", fmt[Max[norms], 4]];
, {k, kVals}];
Print[""];

(* ================================================================ *)
(* 3: At large T, measure error for many consecutive zeros         *)
(* This shows the DISTRIBUTION of errors, not just single points   *)
(* ================================================================ *)

Print["--- 3: Error distribution near zero #10000 (T ~ 9878) ---"];
Print[""];

k0 = 100;
errDist = Table[
  t1 = N[Im[ZetaZero[n]]];
  t2 = N[Im[ZetaZero[n + 1]]];
  tmid = (t1 + t2)/2;
  sExact = n - nSmooth[tmid];
  fk = primeSum[tmid, k0];
  Abs[fk - sExact]
, {n, 9900, 10100}];

Print["k = ", k0, ", zeros #9900 to #10100 (T ~ 9878):"];
Print["  mean |error| = ", fmt[Mean[errDist], 4]];
Print["  std dev      = ", fmt[StandardDeviation[errDist], 4]];
Print["  max |error|  = ", fmt[Max[errDist], 4]];
Print["  fraction < 0.5: ", fmt[N@Count[errDist, x_ /; x < 0.5]/Length[errDist], 3]];
Print["  fraction < 0.3: ", fmt[N@Count[errDist, x_ /; x < 0.3]/Length[errDist], 3]];
Print["  spacing at T=9878: ", fmt[spacing[9878.], 4]];
Print[""];

(* Same for larger T *)
Print["--- 4: Error distribution near zero #100000 (T ~ 74921) ---"];
Print[""];

errDist2 = Table[
  t1 = N[Im[ZetaZero[n]]];
  t2 = N[Im[ZetaZero[n + 1]]];
  tmid = (t1 + t2)/2;
  sExact = n - nSmooth[tmid];
  fk = primeSum[tmid, k0];
  Abs[fk - sExact]
, {n, 99900, 100100}];

Print["k = ", k0, ", zeros #99900 to #100100 (T ~ 74921):"];
Print["  mean |error| = ", fmt[Mean[errDist2], 4]];
Print["  std dev      = ", fmt[StandardDeviation[errDist2], 4]];
Print["  max |error|  = ", fmt[Max[errDist2], 4]];
Print["  fraction < 0.5: ", fmt[N@Count[errDist2, x_ /; x < 0.5]/Length[errDist2], 3]];
Print["  fraction < 0.3: ", fmt[N@Count[errDist2, x_ /; x < 0.3]/Length[errDist2], 3]];
Print["  spacing at T=74921: ", fmt[spacing[74921.], 4]];
Print[""];

(* ================================================================ *)
(* 5: The real question — does the error GROW with T?              *)
(* S(T) itself grows as O(log T / log log T). Does our approx     *)
(* error grow too?                                                  *)
(* ================================================================ *)

Print["--- 5: Does error grow with T? ---"];
Print[""];
Print["Testing with k=100 at exponentially spaced heights:"];
Print[""];

Print[StringJoin[
  StringPadRight["zero#", 9],
  StringPadRight["T", 12],
  StringPadRight["spacing", 10],
  StringPadRight["|S(T)|", 10],
  StringPadRight["|err|", 10],
  StringPadRight["err/spac", 10],
  StringPadRight["logT/loglogT", 14]
]];
Print[StringJoin[Table["-", 75]]];

Do[
  t1 = N[Im[ZetaZero[nz]]];
  t2 = N[Im[ZetaZero[nz + 1]]];
  tmid = (t1 + t2)/2;
  sExact = nz - nSmooth[tmid];
  fk = primeSum[tmid, k0];
  err = Abs[fk - sExact];
  sp = spacing[tmid];

  Print[StringJoin[
    StringPadRight[ToString[nz], 9],
    StringPadRight[fmt[tmid, 5], 12],
    StringPadRight[fmt[sp, 4], 10],
    StringPadRight[fmt[Abs[sExact], 3], 10],
    StringPadRight[fmt[err, 3], 10],
    StringPadRight[fmt[err/sp, 3], 10],
    StringPadRight[fmt[Log[tmid]/Log[Log[tmid]], 3], 14]
  ]];
, {nz, {10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000}}];

Print[""];
Print["If err ~ const (independent of T), then err/spacing grows as log(T)."];
Print["If err ~ 1/log(T), then err/spacing stays constant."];
Print[""];

Print["Done."];
