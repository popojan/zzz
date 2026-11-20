#!/usr/bin/env wolframscript
(* MAX-SNAP-MARGIN DETECTOR experiment.

   Backward step (zeros -> primes) reads off prime powers from the truncated
   psi' signal by thresholding.  Because the answer is an integer, the readout
   is a QUANTIZER: zero error below the lattice resolution snaps away.  This
   script asks whether a TAPERED / MAJORANT window on the zero-sum widens the
   snap margin (the largest zero noise eps still giving a flawless prime
   readout) versus the hard cutoff (Dirichlet) detector that the loop uses.

   Mechanism under test: the hard cutoff produces a Dirichlet (sinc) kernel
   with large Gibbs sidelobes; near the reach edge those deterministic
   sidelobes -- not Gaussian zero-jitter -- drive false positives (spurious
   peak at a composite) and misses (ringing dip at a prime power).  A
   non-negative / tapered kernel (Fejer, Hann, Gaussian) suppresses the
   sidelobes at the cost of a wider main lobe (=> less reach).  We measure
   both faces: snap margin eps*(X) AND the reach cost.

   Method.  Baseline = first nz Odlyzko zeros (doc/refs/odlyzko/zeros1, ~1e-9,
   effectively exact).  Inject i.i.d. Gaussian jitter of size eps (worst case;
   the real loop error is partly common-mode, which cancels more, so this is a
   stress test).  Detector = the project's envelope detector psiP/env with a
   window w(gamma) on the cos-sum, DC-normalized (Total[w]=nz) so EVERY window
   has the same prime-power peak height -- only the sidelobe structure differs,
   so the same tau threshold is a fair comparison.

   Scoring (external; never used by the detector) is loop-faithful: a FAILURE
   at x is a MISSED PRIME (x prime, r<=tau) or a FALSE POSITIVE at a genuine
   composite (x not a prime power, r>tau).  Higher prime powers p^m (m>=2),
   whose peak fades to ~1/m, are DON'T-CARES -- the loop's known-power
   bookkeeping handles them, exactly as in loop.c.  clean reach = largest X
   whose prefix [2,X] is failure-free.  *)

dir = DirectoryName[$InputFileName];
zfile = FileNameJoin[{dir, "..", "refs", "odlyzko", "zeros1"}];

nz   = ToExpression[Environment["MSNZ"]] /. $Failed -> 1200;   (* zeros used *)
Xmax = ToExpression[Environment["MSX"]]  /. $Failed -> 300;    (* scan integers 2..Xmax *)
If[! IntegerQ[nz], nz = 1200]; If[! IntegerQ[Xmax], Xmax = 300];
gam  = N[ToExpression /@ StringTrim /@ Take[ReadList[zfile, String], nz]];

(* ---- detector pieces (project convention) ---------------------------------*)
env[x_, n_] := (20/4) (LogIntegral[n - 4] Log[x])/x;   (* matched-filter peak scale *)
tau = 0.30;                                            (* prime-power threshold *)

xs    = Range[2, Xmax];
logx  = N[Log[xs]];
baseT = 1 + 1/(xs - xs^3);                            (* 1 + 1/(x-x^3) term, vector *)
amp   = 4/Sqrt[N[xs]];                                (* (4/sqrt x) prefactor, vector *)
envv  = env[N[xs], nz];                               (* envelope, vector over x *)
isPrime = PrimeQ /@ xs;                               (* ground truth (scoring only) *)
isPP    = PrimePowerQ /@ xs;                          (* p^m, m>=1 *)

(* ---- windows w(gamma), DC-normalized to Total = nz (equal peak height) -----*)
g1   = gam/Last[gam];                                 (* gamma in [0,1] *)
norm[w_] := w nz/Total[w];
windows = <|
  "Dirichlet" -> norm[ConstantArray[1., nz]],         (* hard cutoff = current detector *)
  "Fejer"     -> norm[1 - g1],                         (* triangular -> non-negative kernel *)
  "Hann"      -> norm[0.5 (1 + Cos[Pi g1])],
  "Gauss"     -> norm[Exp[-(2.0 g1)^2]]
|>;
wkeys = Keys[windows];

(* clean reach for a precomputed r-vector (r = psiW/env over xs): largest X
   with no missed prime and no composite false positive in [2,X]. *)
cleanReach[r_] := Module[{x = Xmax, det},
  det = Thread[r > tau];
  Do[
   If[(isPrime[[i]] && ! det[[i]]) || (! isPP[[i]] && det[[i]]),
      x = xs[[i]] - 1; Break[]],
   {i, Length[xs]}];
  x];

reachFor[gamN_, wk_] := Module[{M = Cos[Outer[Times, logx, gamN]]},
  cleanReach[(baseT - amp (M . windows[wk]))/envv]];

(* ---- sweep ---------------------------------------------------------------- *)
epsGrid = {0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20};
reps    = 20;
SeedRandom[20260621];

(* raw = list of {window, eps, cleanReach}; M shared across windows per draw *)
raw = Flatten[Table[
   Module[{gamN, M},
    gamN = gam + eps RandomVariate[NormalDistribution[], nz];
    M = Cos[Outer[Times, logx, gamN]];
    Table[{wk, eps, cleanReach[(baseT - amp (M . windows[wk]))/envv]}, {wk, wkeys}]],
   {eps, epsGrid}, {rep, reps}], 2];

medReach[wk_, eps_] := Round[Median[N[Cases[raw, {wk, eps, r_} :> r]]]];
epsStar[wk_, X_]    := Module[{ok = Select[epsGrid, medReach[wk, #] >= X &]},
   If[ok === {}, 0., Max[ok]]];

(* ---- report --------------------------------------------------------------- *)
WriteString["stdout", "# MAX-SNAP-MARGIN detector  (nz=" <> ToString[nz] <>
   " zeros, gamma<=" <> ToString[NumberForm[Last[gam], 6]] <> ", scan 2.." <>
   ToString[Xmax] <> ", reps=" <> ToString[reps] <> ")\n"];
WriteString["stdout",
   "# median clean reach (largest flawless X) vs injected zero noise eps\n\n"];
WriteString["stdout", "eps      " <> StringJoin[StringPadLeft[#, 11] & /@ wkeys] <> "\n"];
Do[WriteString["stdout",
   StringPadRight[ToString[NumberForm[eps, {4, 3}]], 9] <>
   StringJoin[StringPadLeft[ToString[medReach[#, eps]], 11] & /@ wkeys] <> "\n"],
   {eps, epsGrid}];

WriteString["stdout", "\n# snap margin  eps*(X) = largest eps holding median reach >= X\n"];
WriteString["stdout", "X        " <> StringJoin[StringPadLeft[#, 11] & /@ wkeys] <> "\n"];
Do[WriteString["stdout",
   StringPadRight[ToString[X], 9] <>
   StringJoin[StringPadLeft[ToString[NumberForm[epsStar[#, X], {4, 3}]], 11] & /@ wkeys] <> "\n"],
   {X, {40, 80, 120, 160, 200, 250}}];

(* reach at eps=0 (resolution limit) and eps=0.036 (loop's native median err) *)
WriteString["stdout", "\n# reach at eps=0 (resolution limit) and eps=0.036 (loop native error)\n"];
WriteString["stdout", "window      reach@0   reach@.036\n"];
Do[Module[{r0, rl},
   r0 = reachFor[gam, wk];
   rl = Round[Median[N[Table[reachFor[gam + 0.036 RandomVariate[NormalDistribution[], nz], wk],
        {reps}]]]];
   WriteString["stdout", StringPadRight[wk, 12] <> StringPadLeft[ToString[r0], 7] <>
      StringPadLeft[ToString[rl], 13] <> "\n"]],
   {wk, wkeys}];

(* ---- TSV artifact --------------------------------------------------------- *)
s = OpenWrite[FileNameJoin[{dir, "maxsnap-detector.tsv"}]];
WriteString[s, "eps\t" <> StringRiffle[wkeys, "\t"] <> "\n"];
Do[WriteString[s, ToString[NumberForm[eps, {4, 3}]] <> "\t" <>
   StringRiffle[ToString[medReach[#, eps]] & /@ wkeys, "\t"] <> "\n"], {eps, epsGrid}];
Close[s];
WriteString["stdout", "\n# wrote " <> FileNameJoin[{dir, "maxsnap-detector.tsv"}] <> "\n"];
