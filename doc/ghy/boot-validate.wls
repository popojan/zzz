(* boot-validate.wls -- compare boot-validate.tsv (from boot-validate.sh)
   against Odlyzko ground truth. Reports per-height mean/max |error| for
   method B and the bootstrap, and the improvement factor.

   Ground truth expected in doc/refs/odlyzko/ (gitignored):
     zeros1 (first 10^5 zeros), zeros3 (10^12), zeros4 (10^21), zeros5 (10^22).
*)

scriptDir = DirectoryName[ExpandFileName[$InputFileName]];
valFile = FileNameJoin[{scriptDir, "boot-validate.tsv"}];
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
loadTable[name_] := parseDec /@ Select[ReadList[FileNameJoin[{odlyzkoDir, name}], String], numLineQ];

zeros1 = loadTable["zeros1"];
z3 = loadTable["zeros3"]; z4 = loadTable["zeros4"]; z5 = loadTable["zeros5"];

(* tag -> truth[offset]: absolute gamma of zero #(base + offset) *)
truth = <|
  "low" -> Function[off, zeros1[[1000 + off]]],
  "mid" -> Function[off, zeros1[[99900 + off]]],
  "e12" -> Function[off, 267653395647 + z3[[off]]],
  "e21" -> Function[off, 144176897509546973000 + z4[[off]]],
  "e22" -> Function[off, 1370919909931995300000 + z5[[off]]]
|>;

rows = Select[StringSplit[#, "\t"] & /@ ReadList[valFile, String], Length[#] == 4 &];

Do[
 Module[{sub, errsB, errsC, offs},
  sub = Select[rows, #[[1]] == tag &];
  If[sub === {}, Continue[]];
  errs[m_] := Module[{rs = Select[sub, #[[3]] == m &]},
    N[parseDec[#[[4]]] - truth[tag][ToExpression[#[[2]]]]] & /@ rs];
  errsB = errs["B"]; errsC = errs["boot"];
  Print[tag, "  (", Length[errsB], " zeros)"];
  Print["  B    mean|err| = ", fmt[Mean[Abs[errsB]]], "   max = ", fmt[Max[Abs[errsB]]]];
  Print["  boot mean|err| = ", fmt[Mean[Abs[errsC]]], "   max = ", fmt[Max[Abs[errsC]]]];
  Print["  improvement    = ", fmt[Mean[Abs[errsB]]/Mean[Abs[errsC]]], "x"];
 ], {tag, {"low", "mid", "e12", "e21", "e22"}}]
