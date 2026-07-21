#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(**)
(*nt[t_]:=t/(2\[Pi]) Log[t/(2\[Pi] E)]+7/8*)
(*wave[q_,t_]:=-((-1)^Floor[-((t Log[q])/\[Pi])]/(4Sqrt[q]))( (-1+Sqrt[q])-#)((1+Sqrt[q])-#)&@*)
(*(Sqrt[ (-1+Sqrt[q])^2+(4 /\[Pi]) Sqrt[q] Log[q]Abs[-(\[Pi]/Log[q])+Mod[-t+\[Pi]/Log[q],(2 \[Pi])/Log[q]]] ])*)
(*ZeroCountApprox[t_,k_]:=nt@t + Sum[wave[Prime@i,t],{i,k}]*)
(**)
(*Plot[{*)
(*ZeroCountApprox[x, 2],*)
(*ZeroCountApprox[x, 100]},{x,10,45},*)
(*GridLines->{Im@ZetaZero@Range@20,Range[10]+1/2},*)
(*PerformanceGoal->"Speed",Frame->True,PlotLabel->"Quadratic Zero Counting Function Approximation",*)
(*PlotLegends->{"k=2","k=100"},FrameLabel->{"imaginary part t","zero count on the critical line"},*)
(*Background->White*)
(*]*)
(**)


(* ::Input:: *)
(**)
(*fix[x_,w_,k_]:=Module[{d},*)
(*prg=0;*)
(*d=First@First@SortBy[Table[{prg+=1;t,Zeta[1/2+I (x+t)]},{t,Range[-w,w,(2w)/k]}],N@Abs@Last@#&];{x+d,w/2,k}]*)
(*nfix[x_,w_,k_,n_]:=First@Nest[fix@@#&,{x,w,k},n]*)
(**)
(*imz[m_,k_,fix_,wp_]:=Module[{x0=((-11+8 m) \[Pi])/(4 ProductLog[(-11+8 m)/(8 E)]),t},*)
(*t=w/.First@FindRoot[nt[w]+Sum[wave[Prime@i,w],{i,k}]==m-1/2,{w,x0},AccuracyGoal->2,PrecisionGoal->2,WorkingPrecision->wp];*)
(*(*Print@"done";*);*)
(*t=If[fix,nfix[t,1/4,20,16],t];*)
(*{N@Zeta[1/2+I t],N[t,30]}*)
(*]*)


(* ::Input:: *)
(*Clear@"x"*)


(* ::Input:: *)
(*prg=0;*)
(*Print@ProgressIndicator[Dynamic[prg],{0,100}];*)


(* ::Input:: *)
(*Timing@imz[10^12+1,2,True,MachinePrecision]*)


(* ::Input:: *)
(*imzplt[m_,k_,wp_]:=Module[{errt=imz[m,k,False,wp],t},*)
(*t=Last@errt;*)
(*ReImPlot[Zeta[1/2+(t+x)I],{x,-Abs@First@errt,Abs@First@errt},PerformanceGoal->"Quality",PlotLabel->t]]*)


(* ::Input:: *)
(*imzplt[1000000,2,MachinePrecision]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Table[Flatten@{2^HoldForm[#]&@m,#,If[m<24,#[[2]]-#[[1]],Null],N@Zeta[1/2+I #[[2]]]}&@{If[m<24,N@Im@ZetaZero[2^m],Null], Last@imz[2^m,2,False,MachinePrecision]},{m,1,36}]//MatrixForm[#,TableHeadings->{None, {"n", "wolfram", "approx", "error","value"}}]&*)


(* ::Input:: *)
(*Plot[nt@x+wave[2,x]+wave[3,x],{x,10,45},GridLines->{Im@ZetaZero@Range@20,Range[10]+1/2}]*)


(* ::Input:: *)
(*halfp[k_,ex_,lo_,hi_]:=Plot[{Mod[nt@x,1]+Sum[wave[Prime@i,x],{i,k}]+wave[ex,x]},{x,lo,hi},GridLines->{((-11+8 #) \[Pi])/(4 ProductLog[(-11+8 #)/(8 E)])&@(Range@20),{1/2}},Epilog->Point[{#,Mod[nt@#,1]+Sum[wave[Prime@i,#],{i,k}]+wave[ex,#]}&/@Im@ZetaZero@Range@20]]*)


(* ::Input:: *)
(*halfp[1,3,0,70]*)


(* ::Input:: *)
(*Plot[{nt@x,nt@x+Sum[pw[Prime@i,x],{i,100}]},{x,90,100},PerformanceGoal->"Speed"]*)


(* ::Input:: *)
(*f*)


(* ::Input:: *)
(*z=Import["~/github/zeta/cmake-build-release/e22.odlyzko.csv"];*)
(*y=Import["~/github/zeta/cmake-build-release/e22.tsv","Data"];*)
(*x=Import["~/github/zeta/cmake-build-release/e22z.tsv","Data"];*)
(*u=Import["~/github/zeta/cmake-build-release/e22_k100.tsv","Data"];*)
(*v=Import["~/github/zeta/cmake-build-release/e22_k1000.tsv","Data"];*)
(*w=Import["~/github/zeta/cmake-build-release/e22_k100000.tsv","Data"];*)


(* ::Input:: *)
(*SmoothHistogram[{*)
(*z[[2;;,1]]-(x[[All,2]]-z[[1,1]]),*)
(*z[[2;;,1]]-(y[[All,2]]-z[[1,1]]),*)
(*z[[2;;,1]]-(u[[All,2]]-z[[1,1]]),*)
(*z[[2;;,1]]-(v[[All,2]]-z[[1,1]]),*)
(*z[[2;;,1]]-(w[[All,2]]-z[[1,1]])*)
(*},PlotLegends->{"k=0","k=2", "k=100","k=1000"(*, "k=100000"*)},Frame->True,*)
(*FrameLabel->{"err (Odlyzko - Zzz)","density"},*)
(*PlotLabel->Row@{Style["Zzz Error Distribution",Bold],"\n\!\(\*SuperscriptBox[\(10\), \(22\)]\) + 1 through \!\(\*SuperscriptBox[\(10\), \(22\)]\) + \!\(\*SuperscriptBox[\(10\), \(4\)]\)"},Background->White]*)


(* ::Input:: *)
(*N[Im@ZetaZero@1000000,10]*)


(* ::Input:: *)
(*fp2[q_,x_]:=Last@p2[q,(Sqrt[q] (Sqrt[(\[Pi]/Log[q]-\[Pi]/(Sqrt[q] Log[q]))^2+(4 \[Pi] x)/(Sqrt[q] Log[q])]-\[Pi]/Log[q]+\[Pi]/(Sqrt[q] Log[q])) Log[q])/(2 \[Pi])]*)
(*fp3[q_,x_]:=Last@p3[q,(-1+Sqrt[q])/(2 Sqrt[q])-(3 \[Pi] (-1+q))/(2^(2/3) q (-108 \[Pi]^3-(54 \[Pi]^3)/q^(3/2)+(162 \[Pi]^3)/Sqrt[q]+216 \[Pi]^2 x Log[q]+Sqrt[(2916 \[Pi]^6 (-1+q)^3)/q^3+(-108 \[Pi]^3-(54 \[Pi]^3)/q^(3/2)+(162 \[Pi]^3)/Sqrt[q]+216 \[Pi]^2 x Log[q])^2])^(1/3))+(-108 \[Pi]^3-(54 \[Pi]^3)/q^(3/2)+(162 \[Pi]^3)/Sqrt[q]+216 \[Pi]^2 x Log[q]+Sqrt[(2916 \[Pi]^6 (-1+q)^3)/q^3+(-108 \[Pi]^3-(54 \[Pi]^3)/q^(3/2)+(162 \[Pi]^3)/Sqrt[q]+216 \[Pi]^2 x Log[q])^2])^(1/3)/(6 2^(1/3) \[Pi])]*)
(*ffp[q_,t_]:=1/4 (Sign@# fp2[q,Abs[#]]&@(Mod[-t+\[Pi]/Log@q,(2\[Pi])/Log@q]-\[Pi]/Log@q)+3Sign@# fp3[q,Abs[#]]&@(Mod[-t+\[Pi]/Log@q,(2\[Pi])/Log@q]-\[Pi]/Log@q))*)
(*ffp2[q_,t_]:=Sign@# fp2[q,Abs[#]]&@(Mod[-t+\[Pi]/Log@q,(2\[Pi])/Log@q]-\[Pi]/Log@q)*)
(*wave[p_,m_]:=1/(\[Pi] m) Exp[-(1/2)m Log[p]]Sin[t m Log[p]]*)
(*sawp[p_,x_]:=Min@{# (Sqrt[p] Log[p])/(p \[Pi]-Sqrt[p] \[Pi]),((-1+Sqrt[p]) Log[p])/((p-1) \[Pi]) (\[Pi]/Log[p]-#)}&@Mod[x,(2\[Pi])/Log[p],-((-Sqrt[p] \[Pi]+p \[Pi])/(2 p Log[p]))]*)
(*p2[p_,t_]:=Module[{a={0,0},b={((-1+Sqrt[p]) \[Pi])/(2 Sqrt[p] Log[p]),(Sqrt[p] Log[p])/(p \[Pi]-Sqrt[p] \[Pi]) ((-1+Sqrt[p]) \[Pi])/(2 Sqrt[p] Log[p])},c={\[Pi]/Log[p],0}},*)
(*(1-t)^2 a+2(1-t)t b  +t^2 c]*)
(*p3[p_,t_]:=Module[{a={0,0},b={((-1+Sqrt[p]) \[Pi])/(2 Sqrt[p] Log[p]),(Sqrt[p] Log[p])/(p \[Pi]-Sqrt[p] \[Pi]) ((-1+Sqrt[p]) \[Pi])/(2 Sqrt[p] Log[p])},c={\[Pi]/Log[p],0}},*)
(*(1-t)^3 a+3(1-t)^2 t b  +3(1-t)t^2 b+t^3 c]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Plot[{(*ffp[Prime@2,t],*)-Sum[wave[Prime@2,i],{i,1,100}],ffp2[Prime@2,t],-sawp[Prime@2,t]},{t,0,16},*)
(*GridLines->{Flatten@{-((-Sqrt[p] \[Pi]+p \[Pi])/(2 p Log[p]))+(2Range[0,3] \[Pi])/Log[p],(-Sqrt[p] \[Pi]+p \[Pi])/(2 p Log[p])+(2Range[0,3] \[Pi])/Log[p]}/.p->3,{}},*)
(*PlotLegends->{"wave sum for p=3 and m = 1 to 100", "quadratic approximation (zzz)", Row@{"tangent saw approximation (limit m" -> \[Infinity],")"}},*)
(*Frame->True,FrameLabel->{"t","counting"},ImageSize->Large,Background->White,PlotLabel->Row@{*)
(*Style["Prime Waves Approximation for p=3",Bold],*)
(*"\n\nwave(p, m):=\!\(\*FractionBox[\(1\), \(\[Pi]\\\ m\)]\)Exp[ -\!\(\*FractionBox[\(1\), \(2\)]\)m Log[p] ] Sin[ t m Log[p] ]\n"}]*)
