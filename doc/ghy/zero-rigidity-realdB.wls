#!/usr/bin/env wolframscript
(* ZERO-RIGIDITY vs the REAL below-threshold B error -- the decisive nail.

   zero-rigidity-resolution.wls showed (Test 2) that the optimal rigidity
   denoiser corrects a SMOOTH synthetic displacement well (33% at coherence
   Lc=3, 92% at Lc=12) and white noise poorly (11%) -- because S(gamma_n) is a
   wiggly log-correlated field, so a coherent displacement is spectrally
   SEPARABLE from it.  But that used a synthetic d independent of the zeros.

   This script generates the ACTUAL method-B error field below threshold
   (kappa = gap*log X < pi, the band-saturation wall) and asks whether the SAME
   oracle rigidity denoiser removes it.  If yes, zero-rigidity genuinely refines
   B beyond the prime-information floor (a real POSITIVE -- the prior-against is
   wrong).  If d_B is secretly in the rigidity signal subspace (uncorrectable),
   the synthetic positive was too easy.  No false positive allowed.

   Honesty: the denoiser is OUT-OF-SAMPLE and ORACLE-trained (target = clean u);
   oracle = "we know the GUE/zeta process covariance," which we do, so this is
   an upper bound on extractable info -- generous to the positive, as required
   for a decisive negative; for a positive it must still be realized via the
   known covariance (next step).  *)

dir   = DirectoryName[$InputFileName];
zfile = FileNameJoin[{dir, "..", "refs", "odlyzko", "zeros1"}];
all   = N[ToExpression /@ StringTrim /@ ReadList[zfile, String]];   (* machine prec: ~10x faster, zeros1 is only ~10 digits anyway *)

n0 = 24001; n1 = 36000;                      (* a clean mid-range block *)
idx = Range[n0, n1];
truez = all[[idx]];

theta[t_] := Im[LogGamma[1/4 + I t/2]] - t/2 Log[Pi];
N0[t_]    := theta[t]/Pi + 1;
uTrue = MapThread[#1 - N0[#2] &, {N[idx], truez}];
uTrue = uTrue - Mean[uTrue];

gapMid = 2 Pi/Log[truez[[Length[truez]/2 // Round]]/(2 Pi)];

(* method-B's estimate of the unfolded fluctuation at each true zero:
   u_B(gamma) = arg P_X(1/2+i gamma)/pi = -(1/pi) sum_{p^m<=X} (1/m) p^{-m/2} sin(m gamma log p).
   This IS the prime-truncation error model of S(gamma) (band-saturation error),
   computed VECTORIZED over the block -- no bisection.  d_B = u_B - S. *)
primesUpTo[X_] := Select[Range[2, X], PrimeQ];
pmList[X_] := Flatten[Reap[Do[Module[{q = p, m = 1},
     While[q <= X, Sow[{p, m}]; m++; q *= p]], {p, primesUpTo[X]}]][[2]], 1];
uBfield[X_, gv_] := Module[{pm = pmList[X], acc},
  acc = Total[(1/(#[[2]] (#[[1]])^(#[[2]]/2)) Sin[#[[2]] Log[#[[1]]] gv]) & /@ pm];
  -acc/Pi];

(* the SAME oracle denoiser as Test 2: estimate clean u_n from a corrupt window
   (incl. centre), weights LS-fit on (corrupt -> clean) out of sample. *)
mFix = 24;
denoiseCorrected[uCorrupt_, uClean_] := Module[
  {Nu = Length[uClean], m = mFix, lo, hi, rows, F, cut, w, te, d},
  d = uCorrupt - uClean;
  lo = m + 1; hi = Nu - m; rows = Range[lo, hi];
  F = Table[uCorrupt[[n - m ;; n + m]], {n, rows}];
  cut = Floor[0.6 Length[rows]];
  w = LeastSquares[F[[;; cut]], uClean[[rows[[;; cut]]]]];
  te = rows[[cut + 1 ;;]];
  {StandardDeviation[d[[te]]], 1 - StandardDeviation[(F[[cut + 1 ;;]] . w) - uClean[[te]]]/StandardDeviation[d[[te]]]}];

WriteString["stdout", "# ZERO-RIGIDITY vs REAL below-threshold B error  (block n=" <>
   ToString[n0] <> ".." <> ToString[n1] <> ", gamma~" <>
   ToString[NumberForm[truez[[Length[truez]/2 // Round]], 6]] <> ", gap=" <>
   ToString[NumberForm[gapMid, 3]] <> ")\n"];
WriteString["stdout", "# denoiser corrected = 1 - ||denoised-truth||/||B-truth|| (1=fully fixed, 0=blind)\n\n"];
WriteString["stdout", "  X   pi(X)  kappa   SD(d_B)/gap   B_err_t   corrected\n"];

Do[Module[{uB, kappa, sdB, cc},
   kappa = gapMid Log[X];
   uB = uBfield[X, truez]; uB = uB - Mean[uB];
   {sdB, cc} = denoiseCorrected[uB, uTrue];
   WriteString["stdout",
     StringPadRight[ToString[X], 4] <> " " <>
     StringPadRight[ToString[Length[primesUpTo[X]]], 6] <> " " <>
     StringPadRight[ToString[NumberForm[kappa, {4, 2}]], 7] <> " " <>
     StringPadRight[ToString[NumberForm[sdB, {6, 4}]], 13] <> " " <>
     StringPadRight[ToString[NumberForm[sdB gapMid, {6, 4}]], 9] <> " " <>
     ToString[NumberForm[cc, {5, 3}]] <> "\n"]],
  {X, {29, 49, 64, 121, 256}}];

WriteString["stdout", "\n# kappa<pi(=3.14) rows are BELOW the band-saturation wall.\n"];
WriteString["stdout", "# If 'corrected' is materially >0 there, the rigidity prior removes part of the\n"];
WriteString["stdout", "# REAL B error that prime-only estimators provably cannot -> decisive POSITIVE.\n"];
