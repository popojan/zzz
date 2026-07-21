#!/usr/bin/env wolframscript
(* Backward reach of the zeros->primes step (Chebyshev psi' explicit formula).
   How far in x can prime powers be read off the first n zeros, when those zeros
   are known only to absolute accuracy eps?  Extends the 2023 exploration
   ~/Documents/Wolfram/2023-10-27_psi_derivative_primepi.nb (kernel dPsiCosLRH +
   envelope (20/4)Li[n-4]log(x)/x) toward the primes<->zeros bootstrap question.
   See doc/notes/mixed-anchor-seeds.md, band-saturation.md.

   Emits TSV (cols: n  eps  x_max  gamma_n) to stdout and bootstrap-backward.tsv.
   Deterministic (SeedRandom). x_max = largest X s.t. the envelope-normalised peak
   perfectly separates prime-powers from non-prime-powers on [2,X]. *)

SeedRandom[1];
maxN = 1200;
allz = N[Im[ZetaZero[Range[maxN]]], 30];

env[x_, n_]   := (20/4) (LogIntegral[n - 4] Log[x])/x;
psiP[x_, zz_] := 1 + 1/(x - x^3) - (4/Sqrt[x]) Total[Cos[zz Log[x]]];

sepReach[n_, eps_, Xtest_] := Module[{zz, r, ppmin = Infinity, npmax = -Infinity, res = Xtest},
  zz = allz[[1 ;; n]];
  If[eps > 0, zz = zz + eps RandomVariate[NormalDistribution[], n]];
  Catch[Do[
     r = N[psiP[x, zz]/env[x, n]];
     If[PrimePowerQ[x], ppmin = Min[ppmin, r], npmax = Max[npmax, r]];
     If[ppmin <= npmax, Throw[res = x - 1]],
     {x, 2, Xtest}]; res];
  res];

ns   = {100, 200, 400, 800, 1200};
epss = {0, 1*^-6, 1*^-4, 1*^-3, 1*^-2, 3*^-2, 1*^-1};
f[x_] := ToString[CForm[N[x]]];

rows = {{"n", "eps", "x_max", "gamma_n"}};
Do[
  xm = Median[Table[sepReach[n, eps, Min[700, Ceiling[2 allz[[n]]]]], {3}]];
  AppendTo[rows, {ToString[n], f[eps], ToString[xm], ToString[Round[allz[[n]]]]}],
  {n, ns}, {eps, epss}];

writeTSV[path_, rs_] := Module[{s = OpenWrite[path]},
  WriteString[s, StringRiffle[StringRiffle[#, "\t"] & /@ rs, "\n"] <> "\n"]; Close[s]];
WriteString["stdout", StringRiffle[#, "\t"] <> "\n"] & /@ rows;
writeTSV[FileNameJoin[{DirectoryName[$InputFileName], "bootstrap-backward.tsv"}], rows];
