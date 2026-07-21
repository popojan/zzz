#!/usr/bin/env wolframscript

(* ================================================================ *)
(* Admissible test function search for the Weil explicit formula    *)
(* Goal: find h(r) whose cosine transform g(u) approximates the    *)
(* heuristic attenuation att(p,T) = 1 - Exp[-Sqrt[T/p]]           *)
(* See: doc/rigor-backlog.md, Phase 1                              *)
(* ================================================================ *)

(* helper: format a number to n decimal places *)
fmt[x_, n_:5] := ToString[AccountingForm[x, n]];

Print[""];
Print["Admissible Test Function Search"];
Print["================================"];
Print[""];

(* ================================================================ *)
(* PART 1: Inverse FT of the heuristic attenuation                 *)
(* ================================================================ *)

Print["--- Part 1: Heuristic attenuation inverse FT ---"];
Print[""];
Print["Claim: Integral Exp[-a e^{-u/2}] e^{iru} du = -2 a^{2ir} Gamma[-2ir]"];
Print[""];

(* Use moderate a for numerical verification *)
aTest = 5.0;
Print["Numerical verification (a = 5):"];
Do[
  num = NIntegrate[
    Exp[-aTest Exp[-u/2]] Exp[I r u],
    {u, -60, 60},
    WorkingPrecision -> 30, MaxRecursion -> 30
  ];
  sym = N[-2 aTest^(2 I r) Gamma[-2 I r], 30];
  Print["  r = ", fmt[r, 4], "  |error| = ", fmt[Abs[num - sym], 3]];
, {r, {0.5, 1.0, 5.0, 10.0, 50.0}}];

Print[""];
Print["Poles of Gamma(-2iz) approaching Im(z) = 1/2:"];
Do[
  val = Abs[N[Gamma[-2 I (1 + I y)]]];
  Print["  y = ", fmt[y, 4], "  |Gamma(-2i(1+iy))| = ", fmt[val, 6]];
, {y, {0, 0.3, 0.45, 0.49, 0.499}}];

Print[""];
Print["Strip of analyticity: |Im(r)| < 1/2 exactly."];
Print["Weil formula needs: |Im(r)| < 1/2 + eps."];
Print["=> HEURISTIC IS NOT ADMISSIBLE (but sits on the boundary)"];
Print[""];

(* ================================================================ *)
(* PART 2: Admissible candidates                                    *)
(* ================================================================ *)

Print["--- Part 2: Admissible candidates ---"];
Print[""];

T0 = 10.0^12;
L = Log[T0];
alpha = Sqrt[T0];

(* Heuristic attenuation *)
attH[u_] := 1 - Exp[-alpha Exp[-u/2]];

(* A: Gaussian step — h(r) = Erfc[(|r|-T)/(sqrt2 Delta)]/2 *)
(* Effective attenuation: Exp[-Delta^2 u^2 / 2] *)
attA[u_, d_] := Exp[-d^2 u^2 / 2];

(* B: Logistic step — h(r) = 1/(1+Exp[(|r|-T)/Delta]) *)
(* g(x) = Pi*Delta*Sin[xT]/Sinh[Pi*Delta*x] *)
(* Effective attenuation: Pi*Delta*u/Sinh[Pi*Delta*u] *)
attB[u_, d_] := If[Abs[u] < 10^-10, 1.0, Pi d u / Sinh[Pi d u]];

(* C: Gaussian with Delta = c/log(T) — adaptive width *)
(* Effective attenuation: Exp[-c^2 u^2 / (2 log(T)^2)] *)
attC[u_, c_] := Exp[-c^2 u^2 / (2 L^2)];

Print["A. Gaussian   att(u) = Exp[-D^2 u^2/2]       strip: infinite"];
Print["B. Logistic   att(u) = pi*D*u/Sinh[pi*D*u]   strip: pi*D > 1/2"];
Print["C. Adaptive   att(u) = Exp[-c^2 u^2/(2L^2)]  strip: infinite"];
Print["   (L = log T = ", fmt[L, 4], ")"];
Print[""];

(* ================================================================ *)
(* PART 3: Parameter optimization                                   *)
(* ================================================================ *)

Print["--- Part 3: Best-fit parameters at T = 10^12 ---"];
Print[""];

uSample = Table[u, {u, Log[2.0], 2.0 L, 0.2}];
target = attH /@ uSample;

(* A: best Delta *)
errAfn[d_?NumericQ] := Mean[(target - (attA[#, d] & /@ uSample))^2];
resA = Quiet[NMinimize[{errAfn[d], d > 0.001}, {d}, Method -> "NelderMead"]];
dA = d /. resA[[2]];
Print["A. Gaussian:  Delta = ", fmt[dA], "  RMS = ", fmt[Sqrt[resA[[1]]], 4]];

(* B: best Delta with admissibility *)
errBfn[d_?NumericQ] := Mean[(target - (attB[#, d] & /@ uSample))^2];
resB = Quiet[NMinimize[{errBfn[d], d > 1/(2 Pi) + 0.001}, {d}, Method -> "NelderMead"]];
dB = d /. resB[[2]];
Print["B. Logistic:  Delta = ", fmt[dB], "  RMS = ", fmt[Sqrt[resB[[1]]], 4],
  "  strip = ", fmt[Pi dB, 4]];

(* C: best c *)
errCfn[c_?NumericQ] := Mean[(target - (attC[#, c] & /@ uSample))^2];
resC = Quiet[NMinimize[{errCfn[c], c > 0.01}, {c}, Method -> "NelderMead"]];
cC = c /. resC[[2]];
Print["C. Adaptive:  c = ", fmt[cC],
  "  (Delta = ", fmt[cC/L], ")",
  "  RMS = ", fmt[Sqrt[resC[[1]]], 4]];
Print[""];

(* ================================================================ *)
(* PART 4: Comparison table                                         *)
(* ================================================================ *)

Print["--- Part 4: Attenuation at selected primes ---"];
Print[""];

row[label_, u_] := StringJoin[
  StringPadRight[label, 10],
  StringPadRight[fmt[u, 4], 9],
  StringPadRight[fmt[attH[u], 5], 11],
  StringPadRight[fmt[attA[u, dA], 5], 11],
  StringPadRight[fmt[attB[u, dB], 5], 11],
  StringPadRight[fmt[attC[u, cC], 5], 11]
];

Print[StringJoin[
  StringPadRight["p", 10], StringPadRight["log p", 9],
  StringPadRight["Heuristic", 11], StringPadRight["Gauss A", 11],
  StringPadRight["Logist B", 11], StringPadRight["Adapt C", 11]
]];
Print[StringJoin[Table["-", 63]]];

testP = {2, 5, 13, 101, 1009, 7919, 10^5, 10^6, 10^8, 10^10, 10^12, 10^14, 10^18};
Do[
  label = If[p < 10^5, ToString[p], "10^" <> ToString[Round[Log10[p]]]];
  Print[row[label, Log[N[p]]]];
, {p, testP}];
Print[""];

(* ================================================================ *)
(* PART 5: Effective weight of first k primes                       *)
(* ================================================================ *)

Print["--- Part 5: Mean attenuation for first k primes ---"];
Print[""];
Print["All used primes have log p << log T = ", fmt[L, 4]];
Print[""];

Do[
  primeU = Table[N@Log[Prime[i]], {i, k}];
  wH = Mean[attH /@ primeU];
  wC = Mean[attC[#, cC] & /@ primeU];
  wB = Mean[attB[#, dB] & /@ primeU];
  Print["  k=", StringPadRight[ToString[k], 6],
    "  Heuristic=", fmt[wH, 4],
    "  Adaptive=", fmt[wC, 4],
    "  Logistic=", fmt[wB, 4]];
, {k, {100, 1000, 10000}}];
Print[""];

(* ================================================================ *)
(* PART 6: Tail analysis                                            *)
(* ================================================================ *)

Print["--- Part 6: Truncation tail ---"];
Print[""];
Print["Heuristic: att ~ sqrt(T)/p for p >> T => tail converges absolutely"];
Print["Admissible: att ~ 1 for used primes => tail only conditional"];
Print[""];
Print["The heuristic buys absolute tail convergence at the cost of"];
Print["non-admissibility. A rigorous approach must bound the conditional"];
Print["tail using partial summation + PNT under RH."];
Print[""];

(* ================================================================ *)
(* PART 7: Gaussian tradeoff table                                  *)
(* ================================================================ *)

Print["--- Part 7: Gaussian (C) tradeoff vs c ---"];
Print[""];
Print["c     att(p=T)   att(p=T^2)  smooth/spacing"];
Print["-------------------------------------------"];
Do[
  delta = c / L;
  spacing = 2 Pi / L;
  Print[fmt[c, 3], "   ",
    fmt[attC[L, c], 5], "     ",
    fmt[attC[2 L, c], 7], "     ",
    fmt[delta/spacing, 4]];
, {c, {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}}];
Print[""];
Print["Need smooth/spacing < 1 => c < 2 Pi ~ 6.28"];

(* Match heuristic at p = T *)
cMatch = c /. FindRoot[attC[L, c] == attH[L], {c, 1}];
Print["Match at p=T: c = ", fmt[cMatch, 4],
  "  att = ", fmt[attH[L], 4],
  "  smooth/spacing = ", fmt[(cMatch/L)/(2 Pi/L), 4]];
Print[""];

(* ================================================================ *)
(* PART 8: Plots                                                    *)
(* ================================================================ *)

Print["--- Part 8: Plots ---"];
Print[""];
Print["The heuristic att = 1 - exp(-v) where v = sqrt(T/p) is scale-free."];
Print["The Gaussian att = exp(-c^2 (log p)^2 / (2 log(T)^2)) is not."];
Print["Plotting in the natural variable v = sqrt(T/p) exposes this."];
Print[""];

(* --- Plot 1: Universal variable v = sqrt(T/p) --- *)
(* Heuristic: 1 - Exp[-v], universal *)
(* Gaussian (C): att = Exp[-c^2 (log T - 2 log v)^2 / (2 log(T)^2)] *)
(*             = Exp[-c^2 (1 - 2 log v / log T)^2 / 2]              *)
attCv[v_, c_, logT_] := Exp[-c^2 (logT - 2 Log[v])^2 / (2 logT^2)];
(* Logistic (B): needs u = log p = log T - 2 log v *)
attBv[v_, d_, logT_] := Module[{u = logT - 2 Log[v]},
  If[Abs[u] < 10^-10, 1.0, Pi d u / Sinh[Pi d u]]
];

p1 = Plot[
  Evaluate[{
    1 - Exp[-v],
    attCv[v, cMatch, L],
    attBv[v, dB, L]
  }],
  {v, 0.01, 30},
  PlotStyle -> {
    {Thick, Black},
    {Thick, Blue, Dashed},
    {Thick, Red, Dashed}
  },
  PlotLegends -> Placed[LineLegend[
    {Directive[Thick, Black], Directive[Thick, Blue, Dashed], Directive[Thick, Red, Dashed]},
    {"Heuristic: 1-exp(-v)  [universal]",
     "Gaussian (c=" <> ToString[Round[cMatch, 0.01]] <> ", T=10^12)",
     "Logistic (T=10^12)"},
    LegendFunction -> "Frame"
  ], {0.65, 0.35}],
  Frame -> True,
  FrameLabel -> {"v = \[Sqrt](T/p)", "Attenuation weight"},
  PlotLabel -> "Attenuation in natural variable v = \[Sqrt](T/p)",
  PlotRange -> {{0, 30}, {-0.05, 1.1}},
  ImageSize -> 600,
  GridLines -> {{1}, {1 - Exp[-1.0]}},
  GridLinesStyle -> Directive[Gray, Dashed],
  Epilog -> {
    Gray,
    Text[Style["v=1 (p=T)", 9], {1.5, -0.03}],
    Text[Style["1-1/e", 9], {28, 0.66}]
  }
];
Export["doc/attenuation-comparison.pdf", p1];
Print["  Saved: doc/attenuation-comparison.pdf"];

(* --- Plot 2: Heuristic is universal, Gaussian is not --- *)
(* Show Gaussian for multiple T values to expose T-dependence *)
p2 = Plot[
  Evaluate[{
    1 - Exp[-v],
    attCv[v, cMatch, Log[10.0^6]],
    attCv[v, cMatch, Log[10.0^12]],
    attCv[v, cMatch, Log[10.0^24]],
    attCv[v, cMatch, Log[10.0^36]]
  }],
  {v, 0.01, 10},
  PlotStyle -> {
    {Thick, Black},
    {Blue, Dashed},
    {Blue},
    {Blue, DotDashed},
    {Blue, Dotted, Thick}
  },
  PlotLegends -> Placed[LineLegend[
    {Directive[Thick, Black],
     Directive[Blue, Dashed], Directive[Blue],
     Directive[Blue, DotDashed], Directive[Blue, Dotted, Thick]},
    {"1-exp(-v) [heuristic, universal]",
     "Gaussian T=10^6", "Gaussian T=10^12",
     "Gaussian T=10^24", "Gaussian T=10^36"},
    LegendFunction -> "Frame"
  ], {0.7, 0.4}],
  Frame -> True,
  FrameLabel -> {"v = \[Sqrt](T/p)", "Attenuation weight"},
  PlotLabel -> "Heuristic is universal; Gaussian depends on T",
  PlotRange -> {{0, 10}, {-0.05, 1.1}},
  ImageSize -> 600,
  GridLines -> {{1}, None},
  GridLinesStyle -> Directive[Gray, Dashed],
  Epilog -> {Gray, Text[Style["p = T", 9], {1, -0.03}]}
];
Export["doc/attenuation-universality.pdf", p2];
Print["  Saved: doc/attenuation-universality.pdf"];

(* --- Plot 3: Zoom on used primes (still in v-space) --- *)
(* For k=1000 at T=10^12: p_max ~ 7919, v = sqrt(10^12/7919) ~ 11237 *)
(* For k=10000: p_max ~ 104729, v ~ 3090 *)
(* All used primes have v >> 1, deep in the att~1 regime *)
Print[""];
Print["  For used primes (k=1000, T=10^12):"];
Print["    p_1000 = 7919,  v = sqrt(T/p) = ", fmt[Sqrt[T0/7919], 1]];
Print["    p_1 = 2,        v = sqrt(T/p) = ", fmt[Sqrt[T0/2], 1]];
Print["    All used primes have v >> 1, deep in att ~ 1 regime."];
Print[""];

(* --- Plot 4: What happens near v=1 (p~T), the transition zone --- *)
p4 = Plot[
  Evaluate[{
    1 - Exp[-v],
    attCv[v, cMatch, Log[10.0^12]],
    attCv[v, 0.5, Log[10.0^12]],
    attCv[v, 2.0, Log[10.0^12]]
  }],
  {v, 0.01, 5},
  PlotStyle -> {
    {Thick, Black},
    {Thick, Blue},
    {Blue, Dashed},
    {Blue, DotDashed}
  },
  PlotLegends -> Placed[LineLegend[
    {Directive[Thick, Black], Directive[Thick, Blue],
     Directive[Blue, Dashed], Directive[Blue, DotDashed]},
    {"Heuristic", "Gaussian c=0.96", "Gaussian c=0.5", "Gaussian c=2.0"},
    LegendFunction -> "Frame"
  ], {0.7, 0.4}],
  Frame -> True,
  FrameLabel -> {"v = \[Sqrt](T/p)", "Attenuation weight"},
  PlotLabel -> "Transition zone near p = T (v = 1)",
  PlotRange -> {{0, 5}, {-0.05, 1.1}},
  ImageSize -> 600,
  GridLines -> {{1}, None},
  GridLinesStyle -> Directive[Gray, Dashed]
];
Export["doc/attenuation-transition.pdf", p4];
Print["  Saved: doc/attenuation-transition.pdf"];
Print[""];

(* ================================================================ *)
(* SUMMARY                                                          *)
(* ================================================================ *)

Print["--- Summary ---"];
Print[""];
Print["1. BOUNDARY: The heuristic's inverse FT involves Gamma(-2ir),"];
Print["   pole at Im(r)=1/2. Exactly on the admissibility boundary."];
Print[""];
Print["2. GAUSSIAN MATCH: Adaptive Gaussian with c ~ ", fmt[cMatch, 3]];
Print["   matches the heuristic at p=T while being always admissible."];
Print[""];
Print["3. USED PRIMES: For p <= p_k (k~10^4), ALL attenuations give ~1."];
Print["   The attenuation shape is irrelevant for the kept terms."];
Print[""];
Print["4. TAIL: The heuristic's edge is absolute tail convergence."];
Print["   Rigorous path: Gaussian att + tail bound via PNT under RH."];
Print[""];
Print["Done."];
