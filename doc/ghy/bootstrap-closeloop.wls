#!/usr/bin/env wolframscript
(* Closed primes<->zeros bootstrap loop (honest: self-computed noisy zeros are
   threaded back into the backward detector).  Seed = finite primes (<=X0) + first
   n0 exact zeros.  Each iteration:
     backward: detect prime powers from current zeros (psi' kernel + envelope, a
               self-contained confidence stop) -> new X_known;
     forward : locate more zeros from primes<=X_known (method B, marching bracket,
               NO truth) -> extend the zero list with realistic errors.
   Tracks (X_known, n_known).  Runs two forward variants:
     B  = plain method B (eps ~ 0.05 floor);
     W  = B error scaled 1/10, MODELLING the measured Weil-fit gain (3-15x).
   Decides whether the Weil lever tips the loop sub-critical -> super-critical.
   See doc/notes/mixed-anchor-seeds.md, band-saturation.md, bootstrap-loop.wls. *)

SeedRandom[2];

(* ---------- shared ---------- *)
allz = N[Im[ZetaZero[Range[1205]]], 30];           (* truth: seed + bracket + W-model *)
gam[n_] := allz[[n]];
env[x_, n_]   := (20/4) (LogIntegral[n - 4] Log[x])/x;
psiP[x_, zz_] := 1 + 1/(x - x^3) - (4/Sqrt[x]) Total[Cos[zz Log[x]]];

(* backward reach = information-theoretic separability limit: largest X s.t. the
   envelope-normalised peak separates prime-powers from the rest on [2,X] (ideal
   detector; a real fixed-threshold detector is a C-side engineering problem).
   Applied to the NOISY self-computed zeros -> honest reach under forward error. *)
backReach[zz_, Xcap_] := Module[{n = Length[zz], r, ppmin = Infinity, npmax = -Infinity, res = Xcap},
  Catch[Do[
     r = N[psiP[x, zz]/env[x, n]];
     If[PrimePowerQ[x], ppmin = Min[ppmin, r], npmax = Max[npmax, r]];
     If[ppmin <= npmax, Throw[res = x - 1]],
     {x, 2, Xcap}]; res];
  res];

(* forward: locate zero #k from primes<=X, marching bracket off prev (no truth) *)
theta[t_] := Im[LogGamma[1/4 + I t/2]] - t/2 Log[Pi];
N0[t_]    := theta[t]/Pi + 1;
ppow[X_]  := Flatten[Table[{p, m}, {p, Prime[Range[PrimePi[X]]]}, {m, 1, Floor[Log[p, X]]}], 1];
FB[t_, pm_] := N0[t] + Im[Total[(1/(#[[2]] Power[#[[1]], #[[2]] (1/2 + I t)])) & /@ pm]]/Pi;
locateMarch[k_, pm_, prev_] := Module[{gap, lo, hi, mid, tgt = k - 1/2, a, b},
  gap = 2 Pi/Log[Max[prev, 10]/(2 Pi)];
  lo = prev + 0.25 gap; hi = prev + 1.9 gap;             (* next zero ~ prev+gap *)
  a = FB[lo, pm] - tgt; b = FB[hi, pm] - tgt;
  If[a b > 0, lo = prev + 0.05 gap; hi = prev + 2.6 gap; a = FB[lo,pm]-tgt; b = FB[hi,pm]-tgt];
  If[a b > 0, Return[prev + gap]];                       (* fallback: predicted *)
  Do[mid = (lo + hi)/2;
     If[(FB[lo, pm] - tgt) (FB[mid, pm] - tgt) <= 0, hi = mid, lo = mid], {40}];
  (lo + hi)/2];

(* ---------- the loop ---------- *)
runLoop[sharpen_, label_] := Module[
  {X = 30, n0 = 60, zhat, niter = 7, Xcap, pm, ntarget, t, traj = {}},
  zhat = gam /@ Range[n0];                               (* seed: exact zeros *)
  Do[
    (* backward *)
    Xcap = Min[700, Floor[1.5 allz[[Length[zhat]]]]];
    Xnew = backReach[zhat, Xcap];
    Xtrue = backReach[gam /@ Range[Length[zhat]], Xcap];  (* audit: exact-zero reach *)
    X = Max[X, Xnew];
    (* forward: extend zeros to ~1.6x using primes<=X, capped at 1200 *)
    pm = ppow[X];
    ntarget = Min[1200, Ceiling[1.6 Length[zhat]]];
    Do[
      t = locateMarch[k, pm, zhat[[k - 1]]];
      zhat = Append[zhat, gam[k] + (t - gam[k]) sharpen],  (* sharpen models Weil *)
      {k, Length[zhat] + 1, ntarget}];
    AppendTo[traj, {it, X, Length[zhat], Round[N[gam[Length[zhat]]]], Xtrue,
       N[Median[Abs[(zhat - gam /@ Range[Length[zhat]])][[-20 ;;]]], 2]}],
    {it, niter}];
  {label, traj}];

f[x_] := ToString[CForm[N[x]]];
report[res_] := Module[{label = res[[1]], traj = res[[2]]},
  WriteString["stdout", "\n# forward = ", label, "\n"];
  WriteString["stdout", "iter\tX_known\tn_known\tgamma_n\tX_exactzeros\teps_fwd(top20)\n"];
  WriteString["stdout", StringRiffle[f /@ #, "\t"] <> "\n"] & /@ traj;];

rB = runLoop[1.0, "B (plain, eps~0.05)"];
rW = runLoop[0.1, "W (Weil-modelled, eps~0.005)"];
report[rB];
report[rW];

(* persist *)
dir = DirectoryName[$InputFileName];
s = OpenWrite[FileNameJoin[{dir, "bootstrap-closeloop.tsv"}]];
WriteString[s, "variant\titer\tX_known\tn_known\tgamma_n\tX_exactzeros\teps_fwd\n"];
Do[WriteString[s, StringRiffle[Prepend[f /@ row, r[[1]]], "\t"] <> "\n"], {r, {rB, rW}}, {row, r[[2]]}];
Close[s];
