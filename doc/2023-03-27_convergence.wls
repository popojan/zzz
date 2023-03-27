#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*pts2[p_]:={{0,0},{(\[Pi]-\[Pi]/Sqrt[p])/(2 Log[p]),1/(2 Sqrt[p])},{\[Pi]/Log[p],0}}*)
(*pts3[p_]:=Module[{a,b,c},{a,b,c}=pts2@p;{a,b,b,c}]*)
(*p2[p_,t_]:=Dot[{(1-t)^2,2(1-t)t ,t^2},p]*)
(*p3[p_,t_]:=Dot[{(1-t)^3,3(1-t)^2 t ,3(1-t)t^2,t^3},p]*)
(*{Last@p2[pts2[p],1/2],Last@p3[pts3[p],1/2]}*)
(*asympSpace[n_]:=-(((-7+8 (-1+n)) \[Pi])/(4 ProductLog[(-7+8 (-1+n))/(8 E)]))+(((-7+8 n) \[Pi])/(4 ProductLog[(-7+8 n)/(8 E)]))(*FullSimplify*)*)


(* ::Input:: *)
(*Sum[1/4 1/Sqrt[p+c i Log[p]Sqrt[p]],{i,1,k}]//FullSimplify(*Prime Gap Harald Cram\[EAcute]r, assumes RH as we do anyway*)*)
(*Sum[3/8 1/Sqrt[p+i 2 Sqrt[p]],{i,1,k}]//FullSimplify(* Prime Twins Worst Case*)*)


(* ::Input:: *)
(*loSum[k_,p_,c_]:=1/4 (HurwitzZeta[1/2,1+Sqrt[p]/(c Log[p])]-HurwitzZeta[1/2,1+k+Sqrt[p]/(c Log[p])])/Sqrt[c Sqrt[p] Log[p]]*)


(* ::Input:: *)
(*hiSum[k_,p_]:=3/8 (HurwitzZeta[1/2,1+Sqrt[p]/2]-HurwitzZeta[1/2,1+k+Sqrt[p]/2])/(Sqrt[2] p^(1/4))*)


(* ::Input:: *)
(*Plot[{loSum[k+1,2,1]-loSum[k,2,1],hiSum[k+1,2]-hiSum[k,2]},{k,1,10},Frame->True,PlotRange->Full]*)


(* ::Input:: *)
(*D[loSum[k,p,c],k]*)


(* ::Input:: *)
(*D[hiSum[k,p],k]*)


(* ::Input:: *)
(*dLoSum[k_,p_,c_]:=HurwitzZeta[3/2,1+k+Sqrt[p]/(c Log[p])]/(8 Sqrt[c Sqrt[p] Log[p]])*)
(*dHiSum[k_,p_]:=(3 HurwitzZeta[3/2,1+k+Sqrt[p]/2])/(16 Sqrt[2] p^(1/4))*)


(* ::Input:: *)
(*Plot[{dLoSum[k,2,1],dHiSum[k,2]},{k,1,1000},Frame->True,PlotRange->Full]*)


(* ::Input:: *)
(*convergencePlot[p_,lo_,hi_]:=Plot[{*)
(*Evaluate@D[loSum[k+1,p,1]-loSum[k,p,1],k],*)
(*Evaluate@D[hiSum[k+1,p]-hiSum[k,p],k]*)
(*},{k,lo,hi},PlotLegends->{"lo","hi"},ScalingFunctions->{"SignedLog","SignedLog"},ImageSize->Large]*)


(* ::Input:: *)
(*convergencePlot[Prime@100000,0,10^4]*)


(* ::Input:: *)
(*Limit[#,{k->\[Infinity],p->\[Infinity]}]&@{*)
(*Evaluate@D[loSum[k+1,p,c]-loSum[k,p,c],k],*)
(*Evaluate@D[hiSum[k+1,p]-hiSum[k,p],k]*)
(*}*)


(* ::Input:: *)
(*b[pi0_,k0_]:={*)
(*loSum[k+1,Prime@pi0,1]-loSum[k,Prime@pi0,1],*)
(*hiSum[k+1,Prime@pi0]-hiSum[k,Prime@pi0]*)
(*}/.k->k0*)


(* ::Input:: *)
(*b[100,0]//FullSimplify//N*)


(* ::Input:: *)
(*nextPrimeContributionUpperBound[p_]:=hiSum[1,p]-hiSum[0,p]*)


(* ::Input:: *)
(*minPrimes[n_,err_]:=Module[{sol,space},*)
(*space=asympSpace[n];*)
(*sol = x/.FindRoot[err space== nextPrimeContributionUpperBound[x],{x,1}];*)
(*{PrimePi@sol,err space//N}]*)


(* ::Input:: *)
(*nextPrimesContribution[p0_,k_]:=Sum[nextPrimeContributionUpperBound[p],{p,NestList[NextPrime,p0,k]}]*)


(* ::Input:: *)
(*nextPrimesContribution[541,10]//N*)
(*nextPrimesContribution[7919,10]//N*)


(* ::Input:: *)
(*contribM=Table[nextPrimesContribution[Prime[2^i],2^k],{i,10},{k,10}]//N;*)
(*contribM//MatrixForm[#,TableHeadings->{Prime[2^HoldForm@#]&/@Range[10],2^HoldForm@#&/@Range[10]}]&*)


(* ::Input:: *)
(*Limit[nextPrimesContribution[p0,k],{p0->\[Infinity],k->\[Infinity]}]*)


(* ::Input:: *)
(*contribM//ListPlot3D[#,AxesLabel->{"initial prime", "next primes count", "next primes contribution"}]&*)


(* ::Input:: *)
(*Plot[{*)
(*nextPrimesContribution[Prime[2^10],Round[2^k]],*)
(*nextPrimesContribution[Prime[Round[2^k]],2^10]*)
(*},{k,1,10},PerformanceGoal->"Speed",GridLines->Automatic]*)


(* ::Input:: *)
(*nextPrimesContributionWorst[p0_,k_]:=Sum[nextPrimeContributionUpperBound[p+2i],{i,k}]*)


(* ::Input:: *)
(*Timing@Limit[nextPrimesContributionWorst[p,k],{k->\[Infinity],p->\[Infinity] }]*)
