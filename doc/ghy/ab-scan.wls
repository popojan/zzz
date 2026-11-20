#!/usr/bin/env wolframscript

(* ================================================================ *)
(* R2': A vs B finder accuracy scan                                  *)
(*                                                                  *)
(* Both methods use only primes. Measure zero-location error vs the  *)
(* Wolfram ZetaZero[] reference across ordinal ranges and k values.  *)
(*                                                                  *)
(* A = zzz default (heuristic damping)                               *)
(* B = zzz -G   (GHY partial Euler P_X with X = p_k)                *)
(*                                                                  *)
(* Report: per (ordinal, k), the signed error  found - reference.    *)
(* ================================================================ *)

fmt[x_, d_:7] := ToString[NumberForm[N[x], {12, d}]];

runZzz[flag_, k_, n_, offset_:0] := Module[{cmd, out},
  cmd = "./build/zzz " <> flag <> " -k " <> ToString[k] <>
        " " <> ToString[n] <> " " <> ToString[offset];
  out = RunProcess[{"bash", "-c", cmd}, "StandardOutput"];
  ToExpression[StringTrim[out]]
];

Print[""];
Print["A vs B finder accuracy scan"];
Print["============================"];
Print[""];

ordinals = {10, 100, 1000, 10000, 100000};
kVals    = {50, 100, 200, 500, 1000};

Print["Fetching Wolfram references..."];
refs = AssociationMap[Function[n, N[Im[ZetaZero[n]], 15]], ordinals];
Do[Print["  #", n, "  ", fmt[refs[n]]], {n, ordinals}];
Print[""];

Print["Running zzz (A) and zzz -G (B) for each (ordinal, k)..."];
results = Association[];
Do[
  Do[
    a = runZzz["",   k, n];
    b = runZzz["-G", k, n];
    eA = a - refs[n];
    eB = b - refs[n];
    results[{n, k}] = <|"A" -> a, "B" -> b, "eA" -> eA, "eB" -> eB|>;
  , {k, kVals}];
, {n, ordinals}];

Print[""];
Print["Results table (signed error, found - reference):"];
Print[""];
header = StringJoin[
  StringPadLeft["ordinal", 10],
  Table[StringPadLeft["k=" <> ToString[k], 14] <> " " <> StringPadLeft["", 14], {k, kVals}]
];
Print["                 ", StringJoin[
  Table[StringPadLeft["k=" <> ToString[k] <> " A/B", 28], {k, kVals}]
]];
Print[StringJoin[Table["-", 10 + 28 * Length[kVals]]]];

Do[
  row = StringPadLeft["#" <> ToString[n], 10];
  Do[
    r = results[{n, k}];
    ea = r["eA"];
    eb = r["eB"];
    row = row <> StringPadLeft[fmt[ea, 5], 14] <> " " <>
                 StringPadLeft[fmt[eb, 5], 13];
  , {k, kVals}];
  Print[row];
, {n, ordinals}];
Print[""];

(* Aggregate: for each k, mean |err| over all ordinals *)
Print["Aggregate |error| averaged over ordinals, per k:"];
Print[""];
Print[StringJoin[
  StringPadRight["k", 8],
  StringPadRight["mean|err_A|", 16],
  StringPadRight["mean|err_B|", 16],
  StringPadRight["max|err_A|", 16],
  StringPadRight["max|err_B|", 16]
]];
Print[StringJoin[Table["-", 72]]];

Do[
  eAs = Abs[results[{#, k}]["eA"]] & /@ ordinals;
  eBs = Abs[results[{#, k}]["eB"]] & /@ ordinals;
  Print[StringJoin[
    StringPadRight[ToString[k], 8],
    StringPadRight[fmt[Mean[eAs]], 16],
    StringPadRight[fmt[Mean[eBs]], 16],
    StringPadRight[fmt[Max[eAs]], 16],
    StringPadRight[fmt[Max[eBs]], 16]
  ]];
, {k, kVals}];
Print[""];

(* ---------------- plots ------------------ *)

(* err vs k at each ordinal *)
absErrA[n_, k_] := Abs[results[{n, k}]["eA"]];
absErrB[n_, k_] := Abs[results[{n, k}]["eB"]];

plotA = ListLogLogPlot[
  Table[Table[{k, absErrA[n, k]}, {k, kVals}], {n, ordinals}],
  Joined -> True, PlotMarkers -> Automatic,
  Frame -> True, FrameLabel -> {"k", "|error|"},
  PlotLabel -> "Heuristic A: finder error vs k",
  PlotLegends -> ("#" <> ToString[#] & /@ ordinals),
  ImageSize -> 600
];

plotB = ListLogLogPlot[
  Table[Table[{k, absErrB[n, k]}, {k, kVals}], {n, ordinals}],
  Joined -> True, PlotMarkers -> Automatic,
  Frame -> True, FrameLabel -> {"k", "|error|"},
  PlotLabel -> "GHY B: finder error vs k",
  PlotLegends -> ("#" <> ToString[#] & /@ ordinals),
  ImageSize -> 600
];

Export["doc/ghy/ab-scan-A.pdf", plotA];
Export["doc/ghy/ab-scan-B.pdf", plotB];
Print["Wrote doc/ghy/ab-scan-A.pdf and ab-scan-B.pdf"];

(* Ratio B/A mean error per k *)
meanA = AssociationMap[Function[k, Mean[Abs[results[{#, k}]["eA"]] & /@ ordinals]], kVals];
meanB = AssociationMap[Function[k, Mean[Abs[results[{#, k}]["eB"]] & /@ ordinals]], kVals];

ratioPlot = ListLogLinearPlot[
  Table[{k, meanB[k] / meanA[k]}, {k, kVals}],
  Joined -> True, PlotMarkers -> Automatic,
  Frame -> True, FrameLabel -> {"k", "mean |err_B| / mean |err_A|"},
  PlotLabel -> "Rigor tax: B / A mean error ratio",
  GridLines -> {None, {1.0}},
  GridLinesStyle -> Directive[Gray, Dashed],
  PlotRange -> {All, {0.5, 3}},
  ImageSize -> 600
];
Export["doc/ghy/ab-scan-ratio.pdf", ratioPlot];
Print["Wrote doc/ghy/ab-scan-ratio.pdf"];

Print["Done."];
