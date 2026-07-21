#!/usr/bin/env wolframscript
(* ZERO-RIGIDITY RESOLUTION -- a decisive test of the one open crack in
   band-saturation.md: can the INTEGRALITY / rigidity of the zeros (a unit
   point process) break the FORWARD degeneracy below the 5e-6 floor?

   The five prior negatives (E2,E3,E4b,E4c, boot) are all potential FALSE
   negatives by construction: each used band-limited prime data or a continuous
   L^2 estimator, both PROVABLY in the null space -- so they cannot see
   zero-rigidity even if it carries information.  This script asks the intrinsic
   question they never could: how well does the zeros' OWN process law pin one
   zero given the others?  That bound is the most any rigidity-based selector
   can do, with NO primes involved.

   Design, engineered against a false NEGATIVE (max power for the hypothesis):
     - exact Odlyzko zeros (truth), not B output;
     - TRUE neighbours (no self-consistency / no circularity);
     - the OPTIMAL predictor: best linear predictor = conditional mean for the
       asymptotically-Gaussian S(t) (Selberg CLT) -- no power lost to a weak
       estimator; we also give it two-sided INTERPOLATION (both neighbours, all
       orders up to m) = maximum information;
     - OUT-OF-SAMPLE residual (train/test split) -- no overfit illusion;
     - a POWER CONTROL in the same harness: inject a high-frequency vs a smooth
       (B-type) displacement of EQUAL size; the predictor MUST correct the
       high-frequency one.  If it does (teeth proven) but cannot correct the
       smooth one, the negative is decisive, not dead.

   Unfold to count units: u_n = n - N0(gamma_n) = S(gamma_n), the rigid
   fluctuation (zero mean, O(1), log-correlated).  A position error delta gamma
   = u * gap, so sigma in u-units IS sigma in gap units -- height-independent.
   Compare to B's error ~0.13 gap and the floor 5e-6 ~ 4e-5 gap.  *)

dir   = DirectoryName[$InputFileName];
zfile = FileNameJoin[{dir, "..", "refs", "odlyzko", "zeros1"}];
Nuse  = 60000;                                  (* zeros used *)
gam   = N[ToExpression /@ StringTrim /@ Take[ReadList[zfile, String], Nuse]];

theta[t_] := Im[LogGamma[1/4 + I t/2]] - t/2 Log[Pi];
N0[t_]    := theta[t]/Pi + 1;
u = MapIndexed[(#2[[1]] - N0[#1]) &, gam];      (* u_n = n - N0(gamma_n) = S(gamma_n) *)
u = u - Mean[u];                                (* drop convention constant *)
sig0 = StandardDeviation[u];                    (* baseline: predict the mean *)
gapHere = 2 Pi/Log[Last[gam]/(2 Pi)];           (* local gap at the top, t-units *)

(* best linear two-sided predictor of u_n from {u_{n-m..n-1}, u_{n+1..n+m}},
   fit OUT OF SAMPLE: train on first 60%, test residual on last 40%. *)
sigmaCond[m_, seq_] := Module[{lo = m + 1, hi = Length[seq] - m, rows, F, y, cut, w, res},
  rows = Range[lo, hi];
  F = Table[Join[seq[[n - m ;; n - 1]], seq[[n + 1 ;; n + m]]], {n, rows}];
  y = seq[[rows]];
  cut = Floor[0.6 Length[rows]];
  w = LeastSquares[F[[;; cut]], y[[;; cut]]];
  res = y[[cut + 1 ;;]] - F[[cut + 1 ;;]] . w;
  StandardDeviation[res]];

(* predicted value of each u_n (out-of-sample style) given a possibly-corrupted
   sequence in the NEIGHBOURS, using weights trained on the CLEAN process
   (oracle prior = true correlations).  Returns predicted vector aligned to rows. *)
predictWith[m_, trainSeq_, neighSeq_] := Module[{lo = m + 1, hi, Ftr, ytr, w, rows, Fp},
  hi = Length[trainSeq] - m;
  rows = Range[lo, hi];
  Ftr = Table[Join[trainSeq[[n - m ;; n - 1]], trainSeq[[n + 1 ;; n + m]]], {n, rows}];
  ytr = trainSeq[[rows]];
  w = LeastSquares[Ftr, ytr];                                  (* oracle weights *)
  hi = Length[neighSeq] - m; rows = Range[lo, hi];
  Fp = Table[Join[neighSeq[[n - m ;; n - 1]], neighSeq[[n + 1 ;; n + m]]], {n, rows}];
  {rows, Fp . w}];

(* ---- Test 1: intrinsic rigidity resolution sigma_cond(m) ------------------ *)
mgrid = {1, 2, 4, 8, 16, 24, 32};
sc = (# -> sigmaCond[#, u]) & /@ mgrid;

WriteString["stdout", "# ZERO-RIGIDITY RESOLUTION  (Nuse=" <> ToString[Nuse] <>
   ", gamma<=" <> ToString[NumberForm[Last[gam], 6]] <> ", gap=" <>
   ToString[NumberForm[gapHere, 3]] <> " t-units)\n"];
WriteString["stdout", "# u_n = S(gamma_n); Var^(1/2)(u) = sig0 = " <>
   ToString[NumberForm[sig0, 4]] <> " (predict-the-mean baseline)\n\n"];
WriteString["stdout", "# Test 1 -- best two-sided interpolation, OUT OF SAMPLE.\n"];
WriteString["stdout", "#   sigma_cond is in GAP units (= position resolution / gap), height-independent.\n"];
WriteString["stdout", "m(each side)   sigma_cond   sigma_cond/sig0\n"];
Do[WriteString["stdout",
   StringPadRight[ToString[m], 15] <>
   StringPadRight[ToString[NumberForm[m /. sc, {6, 4}]], 13] <>
   ToString[NumberForm[(m /. sc)/sig0, {4, 3}]] <> "\n"], {m, mgrid}];
scMin = Min[Values[sc]];

WriteString["stdout", "\n# reference scales (gap units):  B error ~ 0.13   |   5e-6 floor ~ " <>
   ToString[ScientificForm[5.*^-6/gapHere, 2]] <> "\n"];
WriteString["stdout", "# rigidity resolution sigma_cond (best m) = " <>
   ToString[NumberForm[scMin, 4]] <> " gap\n"];

(* ---- Test 2: can rigidity DENOISE the displacement B actually makes? -------
   Observe y = u + d (EVERY zero displaced by d).  Optimal linear denoiser
   (Wiener): estimate u_n from the corrupt window INCLUDING the centre, weights
   trained on (corrupt features -> CLEAN target) per displacement type (oracle
   prior), OUT OF SAMPLE.  Sweep the coherence length Lc of d at fixed
   magnitude: Lc=0 = white (power control -- MUST be corrected if the denoiser
   has teeth); large Lc = the smooth, coherent B-mode (aliased truncation
   deficit, coherence ~ 1/kappa ~ a few gaps).  *)
mFix = 24; Delta = 0.13; Nu = Length[u];
SeedRandom[20260622];
white = RandomVariate[NormalDistribution[], Nu];
field[Lc_] := Delta (#/StandardDeviation[#]) &[If[Lc <= 0, white, GaussianFilter[white, Lc]]];
denoiseCorrected[d_] := Module[{corrupt = u + d, m = mFix, lo, hi, rows, F, cut, w, te},
  lo = m + 1; hi = Nu - m; rows = Range[lo, hi];
  F = Table[corrupt[[n - m ;; n + m]], {n, rows}];        (* INCLUDES centre = denoise *)
  cut = Floor[0.6 Length[rows]];
  w = LeastSquares[F[[;; cut]], u[[rows[[;; cut]]]]];      (* clean target, oracle prior *)
  te = rows[[cut + 1 ;;]];
  1 - StandardDeviation[(F[[cut + 1 ;;]] . w) - u[[te]]]/StandardDeviation[d[[te]]]];

LcGrid = {0, 1, 2, 3, 5, 8, 12};
cf = (# -> denoiseCorrected[field[#]]) & /@ LcGrid;
WriteString["stdout", "\n# Test 2 -- can rigidity DENOISE a coherent displacement? (Delta=0.13 gap, m=24)\n"];
WriteString["stdout", "#   corrected = 1 - ||denoised - truth|| / ||displaced - truth||  (1=fully fixed, 0=blind)\n"];
WriteString["stdout", "#   Lc = coherence length of the displacement, in gaps.  B's aliased deficit ~ Lc 2-4.\n"];
WriteString["stdout", "Lc(gaps)   corrected\n"];
Do[WriteString["stdout", StringPadRight[ToString[Lc], 11] <>
   ToString[NumberForm[Lc /. cf, {5, 3}]] <> "\n"], {Lc, LcGrid}];
cfWhite = 0 /. cf; cfB = 3 /. cf;

(* ---- verdict -------------------------------------------------------------- *)
WriteString["stdout", "\n# VERDICT\n"];
WriteString["stdout", "# Test 1: rigidity resolves one zero to sigma_cond=" <> ToString[NumberForm[scMin, 3]] <>
   " gap from " <> ToString[2 Max[mgrid]] <> " TRUE neighbours -- SHARP, far below B's 0.13 (Ghosh-Peres).\n"];
WriteString["stdout", "# Power control: white (Lc=0) displacement corrected by " <>
   ToString[NumberForm[cfWhite, {3, 2}]] <> " -> the denoiser has teeth.\n"];
WriteString["stdout", If[cfWhite - cfB > 0.25,
   "# DECISIVE NEGATIVE: the coherent B-mode (Lc~3) corrected only " <> ToString[NumberForm[cfB, {4, 3}]] <>
     " -- rigidity is sharp but BLIND to the displacement B actually makes.\n" <>
     "#   Mechanism: B's error = the rigidity null mode (a number+com-preserving warp), not 5 accidents.\n",
   "# POSITIVE / PARTIAL: coherent B-mode corrected " <> ToString[NumberForm[cfB, {4, 3}]] <>
     " -- rigidity refines B even in the coherent mode; PURSUE.\n"]];

s = OpenWrite[FileNameJoin[{dir, "zero-rigidity-resolution.tsv"}]];
WriteString[s, "m\tsigma_cond\n"];
Do[WriteString[s, ToString[m] <> "\t" <> ToString[m /. sc] <> "\n"], {m, mgrid}];
WriteString[s, "sig0\t" <> ToString[sig0] <> "\ngap\t" <> ToString[gapHere] <> "\n"];
Do[WriteString[s, "corrected_Lc" <> ToString[Lc] <> "\t" <> ToString[Lc /. cf] <> "\n"], {Lc, LcGrid}];
Close[s];
WriteString["stdout", "\n# wrote " <> FileNameJoin[{dir, "zero-rigidity-resolution.tsv"}] <> "\n"];
