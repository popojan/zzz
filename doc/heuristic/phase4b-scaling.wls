#!/usr/bin/env wolframscript

(* ================================================================ *)
(* Phase 4b: How does the error depend on T/p_k?                   *)
(* Test at larger T where T >> p_k                                 *)
(* ================================================================ *)

fmt[x_, n_:5] := ToString[AccountingForm[N[x], n]];

Print[""];
Print["Phase 4b: Error scaling with T/p_k"];
Print["===================================="];
Print[""];

nSmooth[t_] := 1 + RiemannSiegelTheta[t]/Pi;

primeSum[t_, k_] := Module[{s = 0.0 + 0.0 I},
  Do[s += Log[1 - Prime[i]^(-1/2 + I t)], {i, k}];
  -Im[s]/Pi
];

(* ================================================================ *)
(* Test 1: Fixed k=100, increasing T                               *)
(* As T grows, T/p_k grows, error should decrease                  *)
(* ================================================================ *)

Print["--- Test 1: Fixed k=100, increasing T ---"];
Print[""];

k0 = 100;
pk = Prime[k0];
Print["k = ", k0, "  p_k = ", pk, "  log(p_k) = ", fmt[Log[N@pk], 3]];
Print[""];

Print[StringJoin[
  StringPadRight["zero#", 8],
  StringPadRight["T", 12],
  StringPadRight["T/p_k", 10],
  StringPadRight["S(T)", 10],
  StringPadRight["F_k(T)", 10],
  StringPadRight["|error|", 10]
]];
Print[StringJoin[Table["-", 60]]];

errByRatio = {};

Do[
  t1 = N[Im[ZetaZero[n]]];
  t2 = N[Im[ZetaZero[n + 1]]];
  tmid = (t1 + t2)/2;
  sExact = n - nSmooth[tmid];
  fk = primeSum[tmid, k0];
  err = Abs[fk - sExact];
  ratio = tmid/pk;
  AppendTo[errByRatio, {ratio, err}];

  If[MemberQ[{10, 50, 100, 200, 500, 1000, 2000, 5000, 10000,
              20000, 50000, 100000}, n],
    Print[StringJoin[
      StringPadRight[ToString[n], 8],
      StringPadRight[fmt[tmid, 6], 12],
      StringPadRight[fmt[ratio, 4], 10],
      StringPadRight[fmt[sExact, 4], 10],
      StringPadRight[fmt[fk, 4], 10],
      StringPadRight[fmt[err, 4], 10]
    ]];
  ];
, {n, Join[Range[1, 200], Range[210, 1000, 10], Range[1050, 5000, 50],
          Range[5100, 10000, 100], Range[10500, 50000, 500],
          Range[51000, 100000, 1000]]}];

Print[""];

(* ================================================================ *)
(* Test 2: Fixed k=100, sweep T — statistics in bins of T/p_k     *)
(* ================================================================ *)

Print["--- Test 2: Error statistics binned by T/p_k ---"];
Print[""];

bins = {{0, 1}, {1, 5}, {5, 10}, {10, 50}, {50, 100}, {100, 500}, {500, Infinity}};
Do[
  {lo, hi} = bin;
  inBin = Select[errByRatio, lo <= #[[1]] < hi &];
  If[Length[inBin] > 0,
    errs = inBin[[All, 2]];
    Print["  T/p_k in [", lo, ", ", hi, "):  n=", Length[inBin],
      "  mean|err|=", fmt[Mean[errs], 4],
      "  max|err|=", fmt[Max[errs], 4],
      "  median=", fmt[Median[errs], 4]];
  ];
, {bin, bins}];
Print[""];

(* ================================================================ *)
(* Test 3: Error at large T with various k                         *)
(* Use T near zero #10000 (~9878) and zero #100000 (~74921)        *)
(* ================================================================ *)

Print["--- Test 3: Large T with various k ---"];
Print[""];

largeTests = {
  {1000, "zero #1000"},
  {5000, "zero #5000"},
  {10000, "zero #10000"},
  {50000, "zero #50000"},
  {100000, "zero #100000"}
};

kVals = {10, 20, 50, 100, 200, 500, 1000};

Do[
  {nz, label} = test;
  t1 = N[Im[ZetaZero[nz]]];
  t2 = N[Im[ZetaZero[nz + 1]]];
  tmid = (t1 + t2)/2;
  sExact = nz - nSmooth[tmid];

  Print[label, "  T ~ ", fmt[tmid, 6], "  S(T) = ", fmt[sExact, 4]];

  Do[
    fk = primeSum[tmid, k];
    err = Abs[fk - sExact];
    Print["  k=", StringPadRight[ToString[k], 5],
      "  T/p_k=", StringPadRight[fmt[tmid/Prime[k], 3], 10],
      "  |error|=", fmt[err, 5]];
  , {k, kVals}];
  Print[""];
, {test, largeTests}];

(* ================================================================ *)
(* Test 4: Error vs T/p_k — scatter plot                           *)
(* ================================================================ *)

Print["--- Test 4: Scatter plot of |error| vs T/p_k ---"];
Print[""];

(* Collect more data points: various k at various T *)
scatterData = {};
Do[
  t1 = N[Im[ZetaZero[nz]]];
  t2 = N[Im[ZetaZero[nz + 1]]];
  tmid = (t1 + t2)/2;
  sExact = nz - nSmooth[tmid];
  Do[
    fk = primeSum[tmid, k];
    err = Abs[fk - sExact];
    AppendTo[scatterData, {tmid/Prime[k], err}];
  , {k, {20, 50, 100, 200, 500}}];
, {nz, {50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000}}];

scatter = ListLogLogPlot[scatterData,
  Frame -> True,
  FrameLabel -> {"T / p_k", "|error|"},
  PlotLabel -> "Truncation error vs T/p_k ratio",
  PlotStyle -> {Blue, PointSize[0.01]},
  ImageSize -> 600,
  GridLines -> {None, {0.5}},
  GridLinesStyle -> Directive[Red, Dashed],
  PlotRange -> All
];
Export["doc/error-vs-ratio.pdf", scatter];
Print["  Saved: doc/error-vs-ratio.pdf"];
Print[""];

(* ================================================================ *)
(* Test 5: Is there a power law?                                   *)
(* ================================================================ *)

Print["--- Test 5: Fitting error ~ (T/p_k)^alpha ---"];
Print[""];

(* Filter to T/p_k > 1 (the regime where the code is designed to work) *)
filtered = Select[scatterData, #[[1]] > 1 &];
logData = Log /@ filtered;
If[Length[logData] > 5,
  fit = LinearModelFit[logData, x, x];
  slope = fit["BestFitParameters"][[2]];
  intercept = fit["BestFitParameters"][[1]];
  Print["  log|error| ~ ", fmt[intercept, 3], " + ", fmt[slope, 3], " log(T/p_k)"];
  Print["  => |error| ~ ", fmt[Exp[intercept], 3], " * (T/p_k)^", fmt[slope, 3]];
  Print["  R^2 = ", fmt[fit["RSquared"], 3]];
];
Print[""];

Print["Done."];
