#!/usr/bin/env wolframscript
(* Verify the primes in a `zzz --loop` checkpoint are CORRECT (every listed
   number is prime) and COMPLETE (no prime up to X_known is missing) -- the
   external, post-hoc check the loop itself deliberately never does.
   Prime entries are the bare-integer lines (headers carry letters, zero
   ordinates carry a '.').

     wolframscript -file doc/ghy/check-loop.wls [zzz-loop.state]                 *)

file = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "zzz-loop.state"];
ls = ReadList[file, String];
ps = ToExpression /@ Select[ls, StringMatchQ[#, DigitCharacter ..] &];
xkl = Select[ls, StringStartsQ[#, "xknown"] &];
xk = If[xkl === {}, Max[ps], ToExpression[StringSplit[First[xkl]][[2]]]];
bound = Max[Max[ps], xk];

wrong = Select[ps, ! PrimeQ[#] &];                       (* listed but composite *)
missing = Complement[Prime /@ Range[PrimePi[bound]], ps]; (* prime <= bound not listed *)
dups = Length[ps] - Length[Union[ps]];
fmt[xs_] := If[xs === {}, "none", Row[{Take[xs, UpTo[10]], If[Length[xs] > 10, " ...", ""]}]];

Print["file         : ", file];
Print["listed primes: ", Length[ps], "  (largest ", Max[ps], ")"];
Print["X_known      : ", xk];
Print["not prime    : ", fmt[wrong]];
Print["missing      : ", fmt[missing], "   (primes <= ", bound, " not listed)"];
Print["duplicates   : ", dups];
ok = wrong === {} && missing === {} && dups == 0;
Print["RESULT       : ", If[ok, "CORRECT and COMPLETE -- exactly the primes up to " <> ToString[bound],
                                 "*** PROBLEM ***"]];
