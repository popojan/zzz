#!/usr/bin/env wolframscript

(* ================================================================ *)
(* R7: A/B/C comparison of counting-function approximations         *)
(*                                                                  *)
(*  A = F_zzz    : current heuristic (main.c, zproxy)               *)
(*  B = F_P      : GHY partial Euler only (zghy)                    *)
(*  C = F_hybrid : GHY partial Euler + local Hadamard (zhybrid)     *)
(*                                                                  *)
(* True counting function between consecutive zeros is the integer  *)
(* index of the preceding zero. Measure |F - N(T)| at T's deep      *)
(* inside each zero interval (midpoint). All three should be close; *)
(* near-zero boundaries, C >> A, B because of the Hadamard factor.  *)
(* ================================================================ *)

binary = "./build/";
zerosFile = "papers/first-20-zeros.txt";
X = 100;

Print[""];
Print["GHY A/B/C comparison"];
Print["====================="];
Print[""];

(* Load zeros *)
zeros = Select[
  ToExpression /@ StringSplit[Import[zerosFile, "String"], "\n"],
  NumericQ
];
Print["Loaded ", Length[zeros], " zeros"];
Print[""];

(* Run a binary and parse TSV output.  cols is the list of column
   indices we care about (1-indexed). *)
runBin[cmd_, cols_] := Module[{out, rows},
  out = RunProcess[{"bash", "-c", cmd}, "StandardOutput"];
  rows = Select[StringSplit[#, "\t"] & /@ StringSplit[out, "\n"], Length[#] >= Max[cols] &];
  Map[#[[cols]] &, ToExpression /@ rows, {2}]
];

(* Evaluate F at midpoints between consecutive zeros
   T_mid[j] = (gamma_j + gamma_{j+1}) / 2
   at these points, N(T) = j *)

midpoints = Table[
  N[(zeros[[j]] + zeros[[j + 1]])/2, 20],
  {j, 1, Length[zeros] - 1}
];

Print["Evaluating 3 methods at ", Length[midpoints], " midpoints..."];
Print[""];

(* For each midpoint T, run three binaries:
     zproxy T T 0 <k>        -- F_zzz (heuristic)
     zghy   T T 0 <X>        -- F_P
     zhybrid <zeros> T T 0 <X> -- F_hybrid
*)

oneBinAt[binary_, args_, col_] := Module[{cmd, out, parts},
  cmd = binary <> " " <> args <> " 2>/dev/null";
  out = RunProcess[{"bash", "-c", cmd}, "StandardOutput"];
  parts = StringSplit[StringTrim[out], "\t"];
  If[Length[parts] >= col, ToExpression[parts[[col]]], Missing[]]
];

evalAll[T_] := Module[{Ts = ToString[DecimalForm[N[T, 15], 15]]},
  {
   oneBinAt["./build/zproxy", Ts <> " " <> Ts <> " 0 " <> ToString[X], 4],
   oneBinAt["./build/zghy",   Ts <> " " <> Ts <> " 0 " <> ToString[X], 4],
   oneBinAt["./build/zhybrid", zerosFile <> " " <> Ts <> " " <> Ts <> " 0 " <> ToString[X], 6]
  }
];

(* Table header *)
Print[StringJoin[
  StringPadRight["#", 5],
  StringPadRight["T_mid", 14],
  StringPadRight["N(T)", 6],
  StringPadRight["F_zzz", 11],
  StringPadRight["F_P", 11],
  StringPadRight["F_hybrid", 11],
  StringPadRight["|err_A|", 10],
  StringPadRight["|err_B|", 10],
  StringPadRight["|err_C|", 10]
]];
Print[StringJoin[Table["-", 93]]];

errsA = {}; errsB = {}; errsC = {};

Do[
  T = midpoints[[j]];
  nTrue = j;
  {fA, fB, fC} = evalAll[T];
  eA = Abs[fA - nTrue];
  eB = Abs[fB - nTrue];
  eC = Abs[fC - nTrue];
  AppendTo[errsA, eA]; AppendTo[errsB, eB]; AppendTo[errsC, eC];
  Print[StringJoin[
    StringPadRight[ToString[j], 5],
    StringPadRight[ToString[NumberForm[N[T], 6]], 14],
    StringPadRight[ToString[nTrue], 6],
    StringPadRight[ToString[NumberForm[fA, 5]], 11],
    StringPadRight[ToString[NumberForm[fB, 5]], 11],
    StringPadRight[ToString[NumberForm[fC, 5]], 11],
    StringPadRight[ToString[NumberForm[eA, 4]], 10],
    StringPadRight[ToString[NumberForm[eB, 4]], 10],
    StringPadRight[ToString[NumberForm[eC, 4]], 10]
  ]];
, {j, 1, Length[midpoints]}];

Print[""];
fmt4[x_] := ToString[NumberForm[N[x], {6, 4}]];
Print["MIDPOINT SUMMARY  (deep inside zero intervals)"];
Print["  method        mean|err|    max|err|"];
Print["  ------        ---------    --------"];
Print["  A heuristic   ", fmt4[Mean[errsA]], "       ", fmt4[Max[errsA]]];
Print["  B P_X only    ", fmt4[Mean[errsB]], "       ", fmt4[Max[errsB]]];
Print["  C hybrid      ", fmt4[Mean[errsC]], "       ", fmt4[Max[errsC]]];
Print[""];

(* Now check just BELOW each zero (stress test for boundary behavior) *)
belowPoints = Table[zeros[[j]] - 0.05, {j, 2, Length[zeros] - 1}];

Print["BELOW-ZERO SUMMARY  (at gamma_j - 0.05, where N(T) = j-1)"];

belowErrsA = {}; belowErrsB = {}; belowErrsC = {};
Do[
  T = belowPoints[[j]];
  nTrue = j;  (* just below gamma_{j+1}, so N=j *)
  {fA, fB, fC} = evalAll[T];
  AppendTo[belowErrsA, Abs[fA - nTrue]];
  AppendTo[belowErrsB, Abs[fB - nTrue]];
  AppendTo[belowErrsC, Abs[fC - nTrue]];
, {j, 1, Length[belowPoints]}];

Print["  method        mean|err|    max|err|"];
Print["  ------        ---------    --------"];
Print["  A heuristic   ", fmt4[Mean[belowErrsA]], "       ", fmt4[Max[belowErrsA]]];
Print["  B P_X only    ", fmt4[Mean[belowErrsB]], "       ", fmt4[Max[belowErrsB]]];
Print["  C hybrid      ", fmt4[Mean[belowErrsC]], "       ", fmt4[Max[belowErrsC]]];
Print[""];

Print["Done."];
