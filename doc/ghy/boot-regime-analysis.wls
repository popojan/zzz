parseDec[s_String] := Module[{t, parts, ip, fp},
  t = StringTrim[s]; parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1, FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  ip + fp];
numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];
z5 = parseDec /@ Select[ReadList["doc/refs/odlyzko/zeros5", String], numLineQ];
base = 1370919909931995300000;
rows = Select[StringSplit[#, "\t"] & /@ ReadList["/tmp/regime.tsv", String], Length[#] == 3 &];
fmt[x_] := ToString[CForm[N[x]]];
Do[Module[{sub, errs},
  sub = Select[rows, #[[2]] == m &];
  errs = N[parseDec[#[[3]]] - (base + z5[[ToExpression[#[[1]]]]])] & /@ sub;
  Print[m, ":  mean|err| = ", fmt[Mean[Abs[errs]]], "   max|err| = ", fmt[Max[Abs[errs]]],
    "   rms = ", fmt[Sqrt[Mean[errs^2]]]];
 ], {m, {"B", "W32R1", "W32R2", "W32R3", "W64R2"}}]
