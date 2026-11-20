#!/usr/bin/env wolframscript
(* Forward accuracy + loop closure for the primes<->zeros bootstrap.
   Forward = method B: F_B(t) = N0(t) + arg P_X(1/2+it)/pi, bisected for zero #n.
   Reads x_max(n) from bootstrap-backward.tsv (run that first), then asks whether
   the primes the BACKWARD step yields (X=x_max(n)) suffice to forward-locate the
   n-th zero to the <=tol the backward step tolerates -- the self-sustaining test.
   See doc/notes/mixed-anchor-seeds.md, band-saturation.md.

   Emits two TSVs next to this script:
     bootstrap-forward.tsv  (cols: n  gamma_n  X  eps_fwd)
     bootstrap-closure.tsv  (cols: n  x_max  eps_fwd_at_xmax  tol  verdict) *)

theta[t_]  := Im[LogGamma[1/4 + I t/2]] - t/2 Log[Pi];
N0[t_]     := theta[t]/Pi + 1;
ppow[X_]   := Flatten[Table[{p, m}, {p, Prime[Range[PrimePi[X]]]}, {m, 1, Floor[Log[p, X]]}], 1];
logPX[t_, pm_] := Total[(1/(#[[2]] Power[#[[1]], #[[2]] (1/2 + I t)])) & /@ pm];
FB[t_, pm_] := N0[t] + Im[logPX[t, pm]]/Pi;

locate[n_, X_, g_] := Module[{pm = ppow[X], lo, hi, mid, gap, tgt = n - 1/2},
  gap = 2 Pi/Log[g/(2 Pi)]; lo = g - 0.9 gap; hi = g + 0.9 gap;
  If[(FB[lo, pm] - tgt) (FB[hi, pm] - tgt) > 0, Return[$Failed]];
  Do[mid = (lo + hi)/2;
     If[(FB[lo, pm] - tgt) (FB[mid, pm] - tgt) <= 0, hi = mid, lo = mid], {44}];
  (lo + hi)/2];

allz = N[Im[ZetaZero[Range[1210]]], 30];
gam[n_] := allz[[n]];
f[x_] := ToString[CForm[N[x]]];
dir = DirectoryName[$InputFileName];
writeTSV[path_, rs_] := Module[{s = OpenWrite[path]},
  WriteString[s, StringRiffle[StringRiffle[#, "\t"] & /@ rs, "\n"] <> "\n"]; Close[s]];

(* x_max(n) at eps=0 from the backward run *)
bw = Rest[Import[FileNameJoin[{dir, "bootstrap-backward.tsv"}], "TSV"]];
toI[c_] := If[StringQ[c], ToExpression[c], c];           (* "100" or 100 -> 100 *)
isZero[c_] := (StringQ[c] && (c == "0." || c == "0")) || (NumberQ[c] && c == 0);
xmax = Association[
   Cases[bw, r_ /; isZero[r[[2]]] :> (toI[r[[1]]] -> toI[r[[3]]])]];
ns = Sort[Keys[xmax]];
Print["x_max(n) read from bootstrap-backward.tsv: ", Normal[xmax]];

(* forward accuracy table over a grid of cutoffs X *)
Xs = {20, 50, 100, 200, 400};
fwd = {{"n", "gamma_n", "X", "eps_fwd"}};
Do[
  t = locate[n, X, gam[n]];
  AppendTo[fwd, {ToString[n], f[gam[n]], ToString[X],
     If[t === $Failed, "fail", f[Abs[t - gam[n]]]]}],
  {n, ns}, {X, Xs}];
writeTSV[FileNameJoin[{dir, "bootstrap-forward.tsv"}], fwd];

(* loop closure: forward eps at X=x_max(n), median over n-1,n,n+1 *)
tol = 0.02;
clo = {{"n", "x_max", "eps_fwd_at_xmax", "tol", "verdict"}};
Do[
  X = xmax[n];
  es = Table[Module[{t = locate[k, X, gam[k]]},
        If[t === $Failed, 99., Abs[t - gam[k]]]], {k, n - 1, n + 1}];
  ef = Median[es];
  v  = Which[ef <= tol, "grows", ef <= 2.5 tol, "marginal", True, "sub-critical"];
  AppendTo[clo, {ToString[n], ToString[X], f[ef], f[tol], v}],
  {n, ns}];
writeTSV[FileNameJoin[{dir, "bootstrap-closure.tsv"}], clo];

WriteString["stdout", "# forward accuracy eps(X,n)\n"];
WriteString["stdout", StringRiffle[#, "\t"] <> "\n"] & /@ fwd;
WriteString["stdout", "\n# loop closure (forward eps at X=x_max(n) vs tol)\n"];
WriteString["stdout", StringRiffle[#, "\t"] <> "\n"] & /@ clo;
