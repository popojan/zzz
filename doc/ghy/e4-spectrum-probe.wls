(* e4-spectrum-probe.wls -- Step 2 of the Weil-formula program: the go/no-go
   singular-spectrum probe (doc/notes/superres-multi-cutoff.md).

   Question: linearized at the TRUE zero configuration, how much of the
   window-zero perturbation space is recoverable from band-limited
   functionals, and at what precision cost?

   The set of available functionals is (at best) the Paley-Wiener space
   PW_a, a = log X, with reproducing kernel K(x,y) = sin(a(x-y))/(pi(x-y)).
   The linearized data for a perturbation {dg_j} is phi'(gamma_j) dg_j, whose
   Riesz representers d/dy K(., gamma_j) have Gram matrix

     G_jk = -s''(gamma_j - gamma_k),  s(z) = sin(a z)/(pi z),
     -s''(z) = (1/pi) [ a^2 sin(az)/z + 2a cos(az)/z^2 - 2 sin(az)/z^3 ],
     G_jj = a^3/(3 pi).

   Eigen-spectrum of G = squared singular values of the ideal design:
   the number of eigenvalues above lambda_1 * 10^(-2P) is the number of
   perturbation directions recoverable with P-digit data. Sqrt of
   diag(G^-1) is the per-zero error amplification (unit-norm functional
   data error -> position error); central zeros vs edge guards separate.

   This is the NECESSARY-condition bound: it ignores far-zero leakage
   (admissibility/window engineering), which only makes things harder.
   If the spectrum says no even here, the program is dead; if yes, the
   sqrt-diag numbers are the engineering targets for step 3.

   Configurations probed (positions = true Odlyzko zeros around 10^22+5000):
     A  a = log p_600   (proxy regime, gap*a ~ 1.20), window +-3   (M ~ 42)
     B  a = log p_600,  window +-6                                  (M ~ 84)
     C  a = log p_10000 (gap*a ~ 1.66), window +-6
     D  positions rescaled to gap*a = 0.91 (the true 10^36, k=10^4 regime),
        a = log p_10000, same M as B
*)

prec = 200;

parseDec[s_String] := Module[{t, parts, ip, fp},
  t = StringTrim[s]; parts = StringSplit[t, "."];
  ip = FromDigits[parts[[1]]];
  fp = If[Length[parts] > 1, FromDigits[parts[[2]]]/10^StringLength[parts[[2]]], 0];
  ip + fp];
numLineQ[s_String] := StringMatchQ[StringTrim[s], RegularExpression["-?[0-9]+\\.[0-9]+"]];
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];
z5 = parseDec /@ Select[ReadList[FileNameJoin[{scriptDir, "..", "refs", "odlyzko", "zeros5"}], String], numLineQ];
Print[Length[z5], " zeros loaded (10^22 block)"];

center = 5000;
T0 = z5[[center]];
relAll = N[#- T0, prec] & /@ z5;   (* offsets, exact differences first *)

(* Gram entry, stable small-z via series *)
gram[a_, z_] := Which[
  z == 0, a^3/(3 Pi),
  Abs[a z] < 1/100, a^3/(3 Pi) - a^5 z^2/(10 Pi) + a^7 z^4/(168 Pi),
  True, (a^2 Sin[a z]/z + 2 a Cos[a z]/z^2 - 2 Sin[a z]/z^3)/Pi];

probe[label_, a_, xs_] := Module[{M = Length[xs], G, evals, l1, S, n40, n80,
    dinv, ctr, gap, need01},
  G = Table[gram[a, xs[[j]] - xs[[k]]], {j, M}, {k, M}];
  evals = Reverse[Sort[Eigenvalues[G]]];
  l1 = evals[[1]];
  gap = N[(xs[[-1]] - xs[[1]])/(M - 1)];
  S = N[2 a (xs[[-1]] - xs[[1]])/2/Pi];
  n40 = Count[evals, e_ /; e > l1 10^-40];
  n80 = Count[evals, e_ /; e > l1 10^-80];
  Print["== ", label, ":  M=", M, "  a=", N[a, 6], "  mean gap=", gap,
    "  gap*a=", N[gap a, 4], "  Shannon S=", NumberForm[S, 4]];
  Print["   log10(lam_n/lam_1) at n=S,S+8,S+16,...:"];
  Print["   ", Table[{n, NumberForm[N[Log10[evals[[n]]/l1]], 5]},
     {n, Select[Join[{Ceiling[S]}, Ceiling[S] + 8 Range[1, 12], {M}], # <= M &]}]];
  Print["   usable ranks: ", n40, "/", M, " above 1e-40,   ", n80, "/", M, " above 1e-80"];
  dinv = Sqrt[Diagonal[Inverse[G]]];
  ctr = Select[Range[M], Abs[xs[[#]] - xs[[(M + 1)/2 // Ceiling]]] < (xs[[-1]] - xs[[1]])/6 &];
  Print["   sqrt diag G^-1 (error amplification, data->position):"];
  Print["     central third: min=", ScientificForm[N[Min[dinv[[ctr]]]], 3],
    "  median=", ScientificForm[N[Median[dinv[[ctr]]]], 3]];
  Print["     edges:         max=", ScientificForm[N[Max[dinv]], 3]];
  Print["   data accuracy for central zero to 0.1 gap: ",
    ScientificForm[N[0.1 gap/Median[dinv[[ctr]]]], 3],
    "   to 0.01 gap: ", ScientificForm[N[0.01 gap/Median[dinv[[ctr]]]], 3]];
  evals];

a600 = N[Log[Prime[600]], prec];     (* 8.3889 *)
a1e4 = N[Log[Prime[10000]], prec];   (* 11.5593 *)

win[xs_, w_] := Select[xs, Abs[#] <= w &];

xsA = win[relAll, 3];
xsB = win[relAll, 6];
probe["A: proxy k=600, +-3", a600, xsA];
probe["B: proxy k=600, +-6", a600, xsB];
probe["C: e22 k=1e4, +-6", a1e4, xsB];

(* D: rescale the true configuration to the 10^36 regime: gap*a = 0.9087
   (gamma_1e36 ~ 8.103e34, mean gap = 2 pi/log(T/2pi) = 0.080066, a = a1e4) *)
gap36 = N[2 Pi/Log[8103080891298886832 10^16/(2 Pi)], prec];
scale = gap36/(N[(xsB[[-1]] - xsB[[1]])/(Length[xsB] - 1), prec]);
xsD = scale # & /@ xsB;
probe["D: 10^36 regime (rescaled), k=1e4", a1e4, xsD];
