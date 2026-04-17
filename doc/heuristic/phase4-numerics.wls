#!/usr/bin/env wolframscript

(* ================================================================ *)
(* Phase 4: Numerical exploration of truncation error               *)
(* Compute S(T) exactly and compare with truncated prime sum        *)
(* ================================================================ *)

fmt[x_, n_:5] := ToString[AccountingForm[N[x], n]];

Print[""];
Print["Phase 4: Truncation Error — Numerical Exploration"];
Print["=================================================="];
Print[""];

(* ================================================================ *)
(* 4.0: Definitions                                                 *)
(* ================================================================ *)

(* Smooth counting function N_0(T) = (T/2pi) log(T/2pi e) + 7/8 *)
n0[t_] := t/(2 Pi) Log[t/(2 Pi E)] + 7/8;

(* Exact smooth part via Riemann-Siegel theta *)
nSmooth[t_] := 1 + RiemannSiegelTheta[t]/Pi;

(* S(T) = N(T) - nSmooth(T), where N(T) counts zeros with 0 < gamma <= T *)
(* Between the n-th and (n+1)-th zero: N(T) = n *)
(* At a zero: N(T) jumps by 1 *)

(* Truncated prime sum (no attenuation, since att ~ 1 for used primes) *)
primeSum[t_, k_] := Module[{s = 0.0 + 0.0 I},
  Do[
    p = Prime[i];
    s += Log[1 - p^(-1/2 + I t)];
  , {i, k}];
  -Im[s]/Pi
];

(* Same but with the actual code's attenuation 1-exp(-2t/sqrt(p)) *)
primeSumAtt[t_, k_] := Module[{s = 0.0 + 0.0 I, att},
  Do[
    p = Prime[i];
    att = 1 - Exp[-2 t / Sqrt[N@p]];
    s += att * Log[1 - p^(-1/2 + I t)];
  , {i, k}];
  -Im[s]/Pi
];

(* ================================================================ *)
(* 4.1: Compute S(T) at known zero locations                       *)
(* ================================================================ *)

Print["--- 4.1: S(T) at sample points ---"];
Print[""];

(* Pick T values between consecutive zeros *)
(* N(T) = n for T between zero #n and zero #(n+1) *)

Print["First 10 zeros for reference:"];
Do[
  z = Im[ZetaZero[n]];
  Print["  zero #", n, " at T = ", fmt[z, 8]];
, {n, 1, 10}];
Print[""];

(* Evaluate S(T) at midpoints between consecutive zeros *)
Print["S(T) at midpoints between zeros:"];
Print[StringJoin[
  StringPadRight["n", 5],
  StringPadRight["T_mid", 14],
  StringPadRight["N_0(T)", 10],
  StringPadRight["nSmooth(T)", 12],
  StringPadRight["S(T)", 10]
]];
Print[StringJoin[Table["-", 51]]];

sValues = {};
Do[
  t1 = Im[ZetaZero[n]];
  t2 = Im[ZetaZero[n + 1]];
  tmid = (t1 + t2)/2;
  nExact = n; (* N(tmid) = n since tmid is between zero n and n+1 *)
  sVal = nExact - nSmooth[tmid];
  AppendTo[sValues, {tmid, sVal}];
  If[n <= 15 || Mod[n, 50] == 0,
    Print[StringJoin[
      StringPadRight[ToString[n], 5],
      StringPadRight[fmt[tmid, 7], 14],
      StringPadRight[fmt[n0[tmid], 5], 10],
      StringPadRight[fmt[nSmooth[tmid], 6], 12],
      StringPadRight[fmt[sVal, 5], 10]
    ]];
  ];
, {n, 1, 200}];
Print["..."];
Print[""];

(* Statistics of S(T) *)
sVals = sValues[[All, 2]];
Print["S(T) statistics over first 200 zero intervals:"];
Print["  mean     = ", fmt[Mean[sVals], 4]];
Print["  std dev  = ", fmt[StandardDeviation[sVals], 4]];
Print["  max |S|  = ", fmt[Max[Abs[sVals]], 4]];
Print["  min S    = ", fmt[Min[sVals], 4]];
Print["  max S    = ", fmt[Max[sVals], 4]];
Print[""];

(* ================================================================ *)
(* 4.2: Prime sum vs S(T) — error as function of k                *)
(* ================================================================ *)

Print["--- 4.2: Truncation error vs number of primes k ---"];
Print[""];

(* Pick a few T values: small, medium, near zero clusters *)
testPoints = {
  {50, (Im[ZetaZero[10]] + Im[ZetaZero[11]])/2, 10},
  {100, (Im[ZetaZero[50]] + Im[ZetaZero[51]])/2, 50},
  {200, (Im[ZetaZero[100]] + Im[ZetaZero[101]])/2, 100},
  {500, (Im[ZetaZero[200]] + Im[ZetaZero[201]])/2, 200}
};

kValues = {5, 10, 20, 50, 100, 200, 500, 1000};

Print["Error |F_k(T) - S(T)| for various k and T:"];
Print[""];

Do[
  {tLabel, tmid, nZero} = tp;
  sExact = nZero - nSmooth[tmid];
  Print["T ~ ", fmt[tmid, 6], " (between zeros #", nZero, " and #", nZero+1,
    ")  S(T) = ", fmt[sExact, 5]];

  Print[StringJoin[
    StringPadRight["  k", 8],
    StringPadRight["F_k(T)", 12],
    StringPadRight["|error|", 12],
    StringPadRight["logT/logP", 12],
    StringPadRight["1/sqrt(P)", 12]
  ]];

  Do[
    fk = primeSum[N[tmid], k];
    err = Abs[fk - sExact];
    pk = Prime[k];
    goldston = Log[tmid]/Log[pk];
    sqrtBound = 1/Sqrt[N@pk];
    Print[StringJoin[
      StringPadRight["  " <> ToString[k], 8],
      StringPadRight[fmt[fk, 5], 12],
      StringPadRight[fmt[err, 5], 12],
      StringPadRight[fmt[goldston, 4], 12],
      StringPadRight[fmt[sqrtBound, 5], 12]
    ]];
  , {k, kValues}];
  Print[""];
, {tp, testPoints}];

(* ================================================================ *)
(* 4.3: Error scaling — systematic sweep                           *)
(* ================================================================ *)

Print["--- 4.3: Error scaling with k ---"];
Print[""];

(* Use T near zero #100, sweep k from 1 to 1000 *)
tmid100 = N[(Im[ZetaZero[100]] + Im[ZetaZero[101]])/2];
s100 = 100 - nSmooth[tmid100];

Print["Fixed T = ", fmt[tmid100, 6], "  S(T) = ", fmt[s100, 5]];
Print[""];

errData = Table[
  fk = primeSum[tmid100, k];
  {k, Abs[fk - s100]}
, {k, Join[Range[1, 50], Range[60, 200, 10], Range[250, 1000, 50]]}];

(* Find where error first drops below various thresholds *)
Do[
  pos = FirstPosition[errData, {_, e_ /; e < thresh}, None];
  If[pos =!= None,
    kMin = errData[[pos[[1]], 1]];
    Print["  |error| < ", thresh, " first at k = ", kMin,
      "  (p_k = ", Prime[kMin], ")"];
  ,
    Print["  |error| < ", thresh, " not achieved with k <= 1000"];
  ];
, {thresh, {1.0, 0.5, 0.2, 0.1, 0.05, 0.01}}];
Print[""];

(* Plot error vs k *)
errPlot = ListLogPlot[errData,
  Joined -> True,
  Frame -> True,
  FrameLabel -> {"k (number of primes)", "|F_k(T) - S(T)|"},
  PlotLabel -> "Truncation error at T ~ " <> ToString[Round[tmid100]],
  PlotRange -> {All, {10^-4, 10}},
  GridLines -> {None, {0.5, 0.1, 0.01}},
  GridLinesStyle -> Directive[Red, Dashed],
  ImageSize -> 600,
  Epilog -> {Red, Dashed,
    Text[Style["0.5", 9, Red], {50, 0.6}],
    Text[Style["0.1", 9, Red], {50, 0.12}],
    Text[Style["0.01", 9, Red], {50, 0.012}]
  }
];
Export["doc/truncation-error-vs-k.pdf", errPlot];
Print["  Saved: doc/truncation-error-vs-k.pdf"];
Print[""];

(* ================================================================ *)
(* 4.4: The exponential sum — magnitude of partial sums            *)
(* ================================================================ *)

Print["--- 4.4: Exponential sum magnitudes ---"];
Print[""];
Print["The key object: |sum_{p <= X} p^{-1/2+iT}|"];
Print[""];

(* For several T, compute cumulative |sum p^{-1/2+iT}| *)
testT = {Im[ZetaZero[100]], 1000.0, 5000.0};

Do[
  Print["T = ", fmt[t, 6], ":"];
  cumSum = 0.0 + 0.0 I;
  Do[
    p = Prime[i];
    cumSum += p^(-1/2 + I t);
    If[MemberQ[{10, 20, 50, 100, 200, 500, 1000}, i],
      Print["  k = ", StringPadRight[ToString[i], 5],
        "  |sum| = ", fmt[Abs[cumSum], 5],
        "  |sum|/sqrt(k) = ", fmt[Abs[cumSum]/Sqrt[i], 4],
        "  sqrt(p_k)/log(p_k) = ", fmt[Sqrt[N@p]/Log[N@p], 4]];
    ];
  , {i, 1, 1000}];
  Print[""];
, {t, testT}];

(* ================================================================ *)
(* 4.5: Does attenuation matter at all?                            *)
(* ================================================================ *)

Print["--- 4.5: Attenuated vs bare sum ---"];
Print[""];

Do[
  {tLabel, tmid, nZero} = tp;
  sExact = nZero - nSmooth[tmid];
  fBare = primeSum[N[tmid], 200];
  fAtt = primeSumAtt[N[tmid], 200];
  Print["T ~ ", fmt[tmid, 5],
    "  S(T) = ", fmt[sExact, 5],
    "  F_bare = ", fmt[fBare, 5],
    "  F_att = ", fmt[fAtt, 5],
    "  |diff| = ", fmt[Abs[fBare - fAtt], 3]];
, {tp, testPoints}];
Print[""];
Print["(Difference should be negligible since att ~ 1 for all p <= p_200)"];
Print[""];

(* ================================================================ *)
(* 4.6: Compare error with Goldston bound                          *)
(* ================================================================ *)

Print["--- 4.6: Actual error vs log(T)/log(P) bound ---"];
Print[""];

(* Compute at multiple T values, fixed k = 100 *)
k0 = 100;
pk0 = Prime[k0];
Print["Fixed k = ", k0, "  (p_k = ", pk0, ")"];
Print[""];
Print[StringJoin[
  StringPadRight["zero#", 7],
  StringPadRight["T", 12],
  StringPadRight["|error|", 12],
  StringPadRight["logT/logP", 12],
  StringPadRight["ratio", 10]
]];
Print[StringJoin[Table["-", 53]]];

ratios = {};
Do[
  tmid = N[(Im[ZetaZero[n]] + Im[ZetaZero[n + 1]])/2];
  sExact = n - nSmooth[tmid];
  fk = primeSum[tmid, k0];
  err = Abs[fk - sExact];
  goldston = Log[tmid]/Log[pk0];
  ratio = If[goldston > 0, err/goldston, 0];
  AppendTo[ratios, ratio];
  If[n <= 10 || Mod[n, 25] == 0,
    Print[StringJoin[
      StringPadRight[ToString[n], 7],
      StringPadRight[fmt[tmid, 6], 12],
      StringPadRight[fmt[err, 5], 12],
      StringPadRight[fmt[goldston, 4], 12],
      StringPadRight[fmt[ratio, 4], 10]
    ]];
  ];
, {n, 1, 200}];

Print[""];
Print["Ratio |error| / (log T / log P):"];
Print["  mean   = ", fmt[Mean[ratios], 4]];
Print["  max    = ", fmt[Max[ratios], 4]];
Print["  median = ", fmt[Median[ratios], 4]];
Print[""];
Print["If ratio << 1, then Goldston bound is loose and there's room"];
Print["for a tighter result."];
Print[""];

Print["Done."];
