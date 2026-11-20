#!/usr/bin/env wolframscript
(* SELF-PAVING primes<->zeros bootstrap -- the "super-critical, demonstrated" run.
   Seed = a finite list of zero ordinates (given data).  The LOOP BODY uses ONLY:
     - the Chebyshev psi' signal  psiP(x) = 1 + 1/(x-x^3) - (4/Sqrt x) sum_k cos(g_k log x)
       evaluated on the loop's OWN (self-computed) zeros, to DETECT prime powers;
     - the smooth counting term N0(t) = theta(t)/pi + 1  (LogGamma only, NO zeta);
     - a partial-Euler sum P_X over the loop's OWN DISCOVERED primes, to LOCATE zeros.
   NO Zeta, NO ZetaZero, NO Prime / PrimeQ / PrimePowerQ / FactorInteger in the loop.
   It paves its own way: zeros -> primes -> more zeros -> more primes ...
   Correctness is checked only AFTER the loop (verification, clearly separated).
   Efficiency is irrelevant; a comfortable kappa margin keeps plain-B accurate enough.
   See doc/notes/mixed-anchor-seeds.md, band-saturation.md, bootstrap-closeloop.wls. *)

(* ===================== SEED DATA (finite, given; not part of the loop) ======= *)
SEEDN = 80;
seedZeros = N[Im[ZetaZero[Range[SEEDN]]], 30];   (* obtained once, externally *)

(* ===================== LOOP PRIMITIVES (zeta-free, primality-free) =========== *)
theta[t_]     := Im[LogGamma[1/4 + I t/2]] - t/2 Log[Pi];
N0[t_]        := theta[t]/Pi + 1;
env[x_, n_]   := (20/4) (LogIntegral[n - 4] Log[x])/x;
psiP[x_, zz_] := 1 + 1/(x - x^3) - (4/Sqrt[x]) Total[Cos[zz Log[x]]];

(* is x a perfect power p^m (m>=2) of an ALREADY-DISCOVERED prime?  (bookkeeping
   against the loop's own prime list -- pure integer arithmetic, not a primality test) *)
isKnownPow[x_, primes_] := Catch[
  Do[q = p p; While[q < x, q *= p]; If[q == x, Throw[True]], {p, primes}]; False];

(* backward: from current zeros, scan x upward and DETECT prime powers from the
   psi' peak (threshold on the envelope).  A new prime = a detected prime-power
   that is not a power of a known prime.  Conservative confidence stop: halt at the
   first ambiguous x that is NOT a known prime power (faded prime OR rising noise) --
   guarantees no false primes corrupt the set.  tau, taulo fixed (no truth). *)
tau = 0.30; taulo = 0.12;
backward[zz_, primes0_, Xcap_] := Module[{n = Length[zz], primes = primes0, r, kp, Xk = 1},
  Catch[Do[
     r  = N[psiP[x, zz]/env[x, n]];
     kp = isKnownPow[x, primes];
     Which[
       r > tau,   (If[! kp && ! MemberQ[primes, x], primes = Append[primes, x]]; Xk = x),
       r > taulo, (If[kp, Xk = x, Throw[Null]]),     (* ambiguous non-power -> stop *)
       True,      Xk = x],                            (* clear composite -> ok *)
     {x, 2, Xcap}]];
  {Sort[primes], Xk}];

(* forward: locate zero #k from primes<=X (their powers), method B, marching bracket
   off the previous self-computed zero and the N0 mean gap.  No truth, no zeta. *)
ppow[primes_, X_] := Flatten[Reap[Do[
     Module[{q = p, m = 1}, While[q <= X, Sow[{p, m}]; m++; q *= p]], {p, Select[primes, # <= X &]}]][[2]], 1];
FB[t_, pm_] := N0[t] + Im[Total[(1/(#[[2]] Power[#[[1]], #[[2]] (1/2 + I t)])) & /@ pm]]/Pi;
locateMarch[k_, pm_, prev_] := Module[{gap, lo, hi, mid, tgt = k - 1/2, a, b},
  gap = 2 Pi/Log[Max[prev, 10]/(2 Pi)];
  lo = prev + 0.25 gap; hi = prev + 1.9 gap; a = FB[lo, pm] - tgt; b = FB[hi, pm] - tgt;
  If[a b > 0, lo = prev + 0.04 gap; hi = prev + 2.7 gap; a = FB[lo, pm] - tgt; b = FB[hi, pm] - tgt];
  If[a b > 0, Return[prev + gap]];
  Do[mid = (lo + hi)/2; If[(FB[lo, pm] - tgt) (FB[mid, pm] - tgt) <= 0, hi = mid, lo = mid], {42}];
  (lo + hi)/2];

(* max zero index resolvable from primes<=X at a comfortable margin kappa>=kmin
   (keeps plain-B accurate enough for the detector; gap=2pi/log(T/2pi)). *)
kmin = 5.0; Nmax = 4500;
nFrontier[X_] := Module[{T = 2 Pi Power[N[X], 2 Pi/kmin]},
  Min[Nmax, Round[(T/(2 Pi)) (Log[T/(2 Pi)] - 1)]]];

(* ===================== THE SELF-PAVING LOOP ================================= *)
zhat = seedZeros; primes = {}; Xknown = 1; niter = 12; traj = {};
Do[
  {primes, Xk} = backward[zhat, primes, Min[400, Ceiling[1.4 Length[zhat]]]];
  Xknown = Max[Xknown, Xk];
  pm = ppow[primes, Xknown];
  nt = Max[Length[zhat], nFrontier[Xknown]];
  Do[zhat = Append[zhat, locateMarch[k, pm, zhat[[k - 1]]]], {k, Length[zhat] + 1, nt}];
  AppendTo[traj, {it, Xknown, Length[primes], Length[zhat]}],
  {it, niter}];

(* ===================== OUTPUT + POST-HOC VERIFICATION ======================= *)
f[x_] := ToString[CForm[N[x]]];
WriteString["stdout", "# self-paving trajectory (loop is zeta-free & primality-free)\n"];
WriteString["stdout", "iter\tX_known\tn_primes\tn_zeros\n"];
WriteString["stdout", StringRiffle[ToString /@ #, "\t"] <> "\n"] & /@ traj;

(* verification (NOT part of the loop): are the discovered primes actually prime,
   and the discovered zeros actually accurate? *)
truePrimesUpto = Select[Range[2, Xknown], PrimeQ];
falsePos = Select[primes, ! PrimeQ[#] &];
missed   = Complement[truePrimesUpto, primes];
nz = Length[zhat];
trueZ = N[Im[ZetaZero[Range[nz]]], 30];
zerr = Abs[zhat - trueZ];
WriteString["stdout", "\n# VERIFICATION (external; not used by the loop)\n"];
WriteString["stdout", "primes_discovered\t" <> ToString[Length[primes]] <> "\n"];
WriteString["stdout", "false_positives\t" <> ToString[falsePos] <> "\n"];
WriteString["stdout", "missed_primes_below_Xknown\t" <> ToString[missed] <> "\n"];
WriteString["stdout", "max_zero_error\t" <> f[Max[zerr]] <> "\n"];
WriteString["stdout", "median_zero_error\t" <> f[Median[zerr]] <> "\n"];
WriteString["stdout", "largest_prime_found\t" <> ToString[Last[primes]] <> "\n"];

dir = DirectoryName[$InputFileName];
s = OpenWrite[FileNameJoin[{dir, "bootstrap-selfpave.tsv"}]];
WriteString[s, "iter\tX_known\tn_primes\tn_zeros\n"];
Do[WriteString[s, StringRiffle[ToString /@ row, "\t"] <> "\n"], {row, traj}];
Close[s];
