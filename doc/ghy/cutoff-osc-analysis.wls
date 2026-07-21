(* cutoff-osc-analysis.wls -- analyse the cutoff sweep produced by cutoff-sweep.sh.
   Tests the parameter-free neighbor-oscillation model for method B's location error:

     Delta(L) = -(1/L) Sum_{j != n} [ Si((gamma_n - gamma_j) L) - (pi/2) sign(gamma_n - gamma_j) ]

   with L = log X = log p_k.  Derivation: F_B(T) = N(T) + (1/pi) Sum_rho Im E_1(i (T - gamma_rho) L),
   Im E_1(iyL) = Si(yL) - (pi/2) sign(y); near the target zero the slope of F_B is L/pi
   (the target's own smeared step), so the n-1/2 crossing shifts by -pi/L times the
   neighbor deficit.  This Delta is exactly what the GHY Z_X factor subtracts, so a high
   model correlation validates the leave-one-out hybrid bootstrap (zzz --boot).

   Findings (2026-06-10 sweep, k = 200..10000):
     n=1000/10000/99950: correlation 0.98, slope ~ 0.9-1.2, model removes 74-80% of rms.
     n=10^12+50:         correlation 0.96, model removes 54% of rms.
     n=100000 degraded (corr 0.64) -- it is the last zero in zeros1, right neighbors missing.
     Plain averaging of gamma over cutoffs buys only 2-4x (nearest-neighbor oscillation
     periods exceed the feasible log-cutoff range); free-frequency LSQ extrapolation is
     ill-conditioned (1 vs 1/L degeneracy).  Direct subtraction (bootstrap) is the winner.

   Ground truth: Odlyzko tables, expected in doc/refs/odlyzko/ (gitignored):
     curl -O https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1   (first 10^5 zeros)
     curl -O https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros3   (10^12 block)
*)

scriptDir = DirectoryName[ExpandFileName[$InputFileName]];
sweepFile = FileNameJoin[{scriptDir, "cutoff-sweep.tsv"}];
If[! FileExistsQ[sweepFile], sweepFile = "/tmp/sweep.tsv"];
odlyzkoDir = FileNameJoin[{scriptDir, "..", "refs", "odlyzko"}];

fmt[x_] := ToString[CForm[N[x]]];

parseDec[s_String] := Module[{t, neg, parts, ip, fp},
  t = StringTrim[s]; neg = StringStartsQ[t, "-"];
  t = StringDelete[t, "-"];
  parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1,
    FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  If[neg, -(ip + fp), ip + fp]];

numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];

sweepRaw = Select[StringSplit[#, "\t"] & /@ ReadList[sweepFile, String], Length[#] == 3 &];
zeros1 = ReadList[FileNameJoin[{odlyzkoDir, "zeros1"}], Number];
z3 = parseDec /@ Select[ReadList[FileNameJoin[{odlyzkoDir, "zeros3"}], String], numLineQ];
base12 = 267653395647;
Print["sweep rows: ", Length[sweepRaw], "   zeros1: ", Length[zeros1], "   zeros3: ", Length[z3]];

(* tag -> {targetGamma, neighborGammas (list incl target), targetIndexInList} *)
mk1[n_, w_] := Module[{lo = Max[1, n - w], hi = Min[Length[zeros1], n + w]},
  {zeros1[[n]], zeros1[[lo ;; hi]], n - lo + 1}];
specs = <|
  "n1000"   -> mk1[1000, 50],
  "n10000"  -> mk1[10000, 50],
  "n99950"  -> mk1[99950, 50],
  "n100000" -> mk1[100000, 50],
  "n1e12"   -> {base12 + z3[[50]], base12 + # & /@ z3[[1 ;; 100]], 50}
|>;

predDelta[L_, tg_, nbrs_, ti_] := Module[{ds},
  ds = N[Table[tg - nbrs[[j]], {j, Delete[Range[Length[nbrs]], ti]}]];
  -(1/L) Total[(SinIntegral[# L] - (Pi/2) Sign[#]) & /@ ds]];

Do[
 Module[{rows, tg, nbrs, ti, ks, Ls, errs, preds, topHalf, errLast, errAvg,
   corr, slope, rmsErr, rmsResid, design, sol, ginf, errFit, gs, freqs, js},
  {tg, nbrs, ti} = specs[tag];
  rows = SortBy[Select[sweepRaw, #[[1]] == tag &], ToExpression[#[[2]]] &];
  ks = ToExpression[#[[2]]] & /@ rows;
  gs = parseDec[#[[3]]] & /@ rows;
  Ls = N[Log[Prime[#]]] & /@ ks;
  errs = N[# - tg] & /@ gs;
  preds = predDelta[#, tg, nbrs, ti] & /@ Ls;

  errLast = Last[errs];
  topHalf = Range[Ceiling[Length[ks]/2], Length[ks]];
  errAvg = N[Mean[gs[[topHalf]]] - tg];

  corr = Correlation[errs, preds];
  slope = errs . preds/(preds . preds);
  rmsErr = Sqrt[Mean[errs^2]];
  rmsResid = Sqrt[Mean[(errs - preds)^2]];

  (* fixed-frequency LSQ: gamma_inf + c/L + sum_{|j|<=3} (a cos + b sin)(d_j L)/L.
     Kept for the record: ill-conditioned (1 vs 1/L degeneracy over feasible L range). *)
  js = Select[DeleteCases[Range[-3, 3], 0], 1 <= ti + # <= Length[nbrs] &];
  freqs = N[tg - nbrs[[ti + #]]] & /@ js;
  design = Table[Join[{1, 1/L}, Flatten[{Cos[# L]/L, Sin[# L]/L} & /@ freqs]], {L, Ls}];
  sol = LeastSquares[design, errs];
  errFit = sol[[1]];

  Print[tag, ":  gamma_true = ", N[tg, 17]];
  Print["  err@kmax         = ", fmt[errLast]];
  Print["  err of avg(topK) = ", fmt[errAvg]];
  Print["  rms(err)         = ", fmt[rmsErr]];
  Print["  model corr       = ", fmt[corr], "   slope(err~pred) = ", fmt[slope]];
  Print["  rms(err - pred)  = ", fmt[rmsResid],
        "   (model removes ", fmt[100 (1 - rmsResid/rmsErr)], "% of rms)"];
  Print["  err of LSQ gamma_inf = ", fmt[errFit], "   (ill-conditioned, see header)"];
 ], {tag, {"n1000", "n10000", "n99950", "n100000", "n1e12"}}]
