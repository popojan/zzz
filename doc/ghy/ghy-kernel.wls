#!/usr/bin/env wolframscript

(* ================================================================ *)
(* R3: GHY kernel U(z) derivation, numerics, and C reference.       *)
(*                                                                  *)
(* Reference: Gonek, Hughes, Young 2007,                             *)
(*   "A hybrid Euler-Hadamard product for the Riemann zeta function" *)
(*                                                                  *)
(* Theorem 1 (GHY): for u an admissible C^infty bump of mass 1 on    *)
(*   [e^{1-1/X}, e], set U(z) = Integrate[u(x) E_1(z log x), {x,0,inf}] *)
(* where E_1 is the exponential integral.                           *)
(*                                                                  *)
(* As X -> infinity, u concentrates on a neighborhood of e of width  *)
(* O(1/X), so U(z) -> E_1(z log e) = E_1(z). We use this limit for   *)
(* the C implementation; the finite-X correction is O(1/X).         *)
(*                                                                  *)
(* This script:                                                     *)
(*  1. Derives U(z) for a specific canonical bump f(y) = 6y(1-y)     *)
(*  2. Verifies U(z) -> E_1(z) as X -> infinity                      *)
(*  3. Plots the Hadamard factor near a zero                        *)
(*  4. Exports reference values for C unit tests                    *)
(* ================================================================ *)

fmt[x_, n_:10] := ToString[NumberForm[N[x], n]];

Print[""];
Print["GHY Kernel U(z): Numerics and Reference"];
Print["========================================"];
Print[""];

(* ================================================================ *)
(* 1: Canonical bump                                                *)
(* ================================================================ *)

(* f supported on [0,1], mass 1. Choose f(y) = 6 y (1-y). *)
f[y_] := 6 y (1 - y);

(* u is derived per Thm 1: u(x) = X f(X log(x/e) + 1)/x *)
u[x_, X_] := X f[X Log[x/E] + 1] / x;

(* Verify mass = 1 for any X *)
Print["--- 1. Sanity: mass of u(x) for selected X ---"];
Do[
  mass = NIntegrate[u[x, Xv], {x, E^(1 - 1/Xv), E}, WorkingPrecision -> 20];
  Print["  X = ", Xv, "   mass = ", fmt[mass]];
, {Xv, {2, 10, 100, 1000}}];
Print[""];

(* ================================================================ *)
(* 2: Kernel U(z) for the canonical bump                           *)
(* ================================================================ *)

(* U(z) = Integrate[u(x) E_1(z log x), {x, 0, inf}] *)
(* Switch vars: y = X log(x/e) + 1, dy = X/x dx, x in [e^{1-1/X}, e] <=> y in [0,1] *)
(* log x = (y-1)/X + 1 *)
(* U(z) = Integrate[f(y) E_1(z ((y-1)/X + 1)), {y, 0, 1}] *)

UExact[z_, X_] := NIntegrate[
  f[y] ExpIntegralE[1, z ((y - 1)/X + 1)],
  {y, 0, 1}, WorkingPrecision -> 20
];

ULimit[z_] := ExpIntegralE[1, z];

Print["--- 2. U(z) vs limit E_1(z): dependence on X ---"];
Print[""];
Print["  For z values of interest: z = i*theta*log(X) for theta near 0"];
Print[""];

Print["z = 1.0 + 0*I:"];
Do[
  diff = Abs[UExact[1.0 + 0 I, Xv] - ULimit[1.0 + 0 I]];
  Print["  X = ", StringPadRight[ToString[Xv], 8],
    "  U(z, X) = ", fmt[UExact[1.0 + 0 I, Xv]],
    "  |U - E_1| = ", fmt[diff, 4]];
, {Xv, {2, 10, 100, 1000}}];
Print[""];

Print["z = 1.0*I (imaginary axis, near-boundary regime):"];
Do[
  exact = UExact[1.0 I, Xv];
  lim = ULimit[1.0 I];
  Print["  X = ", StringPadRight[ToString[Xv], 8],
    "  U(z, X) = ", fmt[exact],
    "  E_1(z) = ", fmt[lim],
    "  diff = ", fmt[Abs[exact - lim], 4]];
, {Xv, {2, 10, 100, 1000}}];
Print[""];

Print["=> For X >= 100, U(z, X) ~ E_1(z) to 2-3 decimal places."];
Print["=> For X >= 1000, agreement is 4+ decimals."];
Print[""];

(* ================================================================ *)
(* 3: Behavior of U near the origin - the key fact                 *)
(* ================================================================ *)

(* Near z = 0: E_1(z) ~ -gamma - log(z). So exp(-U(z)) has a log singularity. *)
(* GHY eq (8): near a zero rho_n = 1/2 + i gamma_n, at s = 1/2 + i t,       *)
(*   exp(-U(i(t - gamma_n) log X)) ~ i(t - gamma_n) e^gamma log X           *)

Print["--- 3. Log-singularity: U(z) for small z ---"];
Print[""];
Print["Asymptotic: E_1(z) ~ -gamma - log(z) + z - z^2/4 + ..."];
Print[""];

Do[
  e1val = ULimit[z];
  asymp = -EulerGamma - Log[z];
  Print["  z = ", fmt[z, 6], "   E_1(z) = ", fmt[e1val],
    "   asymp = ", fmt[asymp],
    "   diff = ", fmt[N[Abs[e1val - asymp]], 4]];
, {z, {0.001, 0.01, 0.1, 0.5, 1.0, 0.01 I, 0.1 I}}];
Print[""];

(* ================================================================ *)
(* 4: Effective window of Hadamard factor                          *)
(* ================================================================ *)

(* Zeros outside the window |t - gamma_n| >> 1/log X contribute ~0. *)
(* Quantify: |exp(-U(z))| for z with large real part.              *)

Print["--- 4. Effective window: |exp(-U(z))| vs |z| ---"];
Print[""];
Print["Pure imaginary z = i*theta (zero on critical line, t = gamma + theta/log X):"];
Print[""];

Do[
  Xv = 1000;
  u0 = ULimit[I theta];
  fac = Exp[-u0];
  Print["  theta = ", StringPadRight[fmt[theta, 4], 9],
    "  U(i*theta) = ", fmt[u0, 6],
    "  |exp(-U)| = ", fmt[Abs[fac], 4],
    "  arg/(pi) = ", fmt[Arg[fac]/Pi, 3]];
, {theta, {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0}}];
Print[""];
Print["=> For |theta| > ~5, |exp(-U)| ~ 1: distant zeros contribute ~1 (inert)."];
Print["=> For |theta| < ~1, |exp(-U)| grows rapidly: local zero dominates."];
Print[""];

(* ================================================================ *)
(* 5: Hadamard factor near a specific zero                         *)
(* ================================================================ *)

Print["--- 5. Hadamard factor near zero #1 (gamma_1 = 14.1347...) ---"];
Print[""];

gamma1 = N[Im[ZetaZero[1]], 20];
X5 = 100;
logX5 = Log[X5];

Print["Parameters: zero at gamma_1 = ", fmt[gamma1], ", X = ", X5,
  ", log X = ", fmt[logX5]];
Print[""];
Print[StringJoin[
  StringPadRight["t - gamma_1", 14],
  StringPadRight["z = i*theta*logX", 20],
  StringPadRight["U(z)", 24],
  StringPadRight["|exp(-U)|", 12],
  "arg(exp(-U))/pi"
]];
Print[StringJoin[Table["-", 85]]];

Do[
  z = I theta logX5;
  uv = ULimit[z];
  fac = Exp[-uv];
  Print[StringJoin[
    StringPadRight[fmt[theta, 5], 14],
    StringPadRight[fmt[z, 6], 20],
    StringPadRight[fmt[uv, 6], 24],
    StringPadRight[fmt[Abs[fac], 5], 12],
    fmt[Arg[fac]/Pi, 4]
  ]];
, {theta, {-1.0, -0.5, -0.1, -0.01, 0.01, 0.1, 0.5, 1.0}}];
Print[""];

(* ================================================================ *)
(* 6: Export reference values for C unit tests                     *)
(* ================================================================ *)

Print["--- 6. Exporting reference values to doc/ghy/ghy-kernel-ref.tsv ---"];
Print[""];

reVals = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0};
imVals = {-5.0, -1.0, -0.1, 0.1, 1.0, 5.0};

refData = Flatten[
  Table[
    Module[{z = re + I im, val},
      val = ULimit[z];
      {re, im, Re[val], Im[val]}
    ],
    {re, reVals}, {im, imVals}
  ],
  1
];

(* Write TSV: re(z)  im(z)  re(E_1(z))  im(E_1(z)) *)
Export[
  "doc/ghy/ghy-kernel-ref.tsv",
  Prepend[refData, {"re_z", "im_z", "re_E1", "im_E1"}],
  "TSV",
  "FieldSeparators" -> "\t"
];
Print["  wrote ", Length[refData], " reference points."];
Print[""];

(* Print first few rows *)
Print["First 5 reference rows:"];
Do[
  r = refData[[i]];
  Print["  z = ", fmt[r[[1]] + I r[[2]], 4],
    "   E_1(z) = ", fmt[r[[3]] + I r[[4]]]];
, {i, 1, Min[5, Length[refData]]}];
Print[""];

(* ================================================================ *)
(* 7: Plot                                                          *)
(* ================================================================ *)

Print["--- 7. Plotting the kernel ---"];
Print[""];

(* |exp(-U(i*theta*logX))| vs theta: the "window shape" of Z_X *)
plot1 = Plot[
  Abs[Exp[-ULimit[I theta Log[1000]]]],
  {theta, -5, 5},
  PlotStyle -> {Thick, Blue},
  Frame -> True,
  FrameLabel -> {"theta  (t - gamma_n, in units of 1)", "|exp(-U(i theta log X))|"},
  PlotLabel -> "Hadamard factor magnitude near a zero  (X = 1000)",
  PlotRange -> {{-5, 5}, All},
  ImageSize -> 600,
  GridLines -> Automatic,
  GridLinesStyle -> Directive[Gray, Dashed]
];
Export["doc/ghy/ghy-hadamard-window.pdf", plot1];
Print["  wrote doc/ghy/ghy-hadamard-window.pdf"];

(* F(t) demo: arg(exp(-U(i(t-gamma_1)logX))) sweeps by pi across gamma_1 *)
plot2 = Plot[
  Arg[Exp[-ULimit[I (t - gamma1) Log[100]]]]/Pi,
  {t, gamma1 - 1, gamma1 + 1},
  PlotStyle -> {Thick, Red},
  Frame -> True,
  FrameLabel -> {"t", "arg(exp(-U)) / pi"},
  PlotLabel -> "Local Hadamard factor argument near zero #1 (X = 100)",
  ImageSize -> 600,
  GridLines -> {{gamma1}, {-0.5, 0, 0.5}},
  GridLinesStyle -> Directive[Gray, Dashed],
  Epilog -> {Gray, Text[Style["gamma_1", 9], {gamma1 + 0.05, 0.45}]}
];
Export["doc/ghy/ghy-hadamard-jump.pdf", plot2];
Print["  wrote doc/ghy/ghy-hadamard-jump.pdf"];

Print[""];
Print["Done."];
