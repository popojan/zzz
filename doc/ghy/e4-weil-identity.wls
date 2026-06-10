(* e4-weil-identity.wls -- Step 1 of the Weil-formula program
   (doc/notes/superres-multi-cutoff.md, "the one principled escape").

   Validates the Riemann-Weil explicit formula as an exact identity against
   Odlyzko's 1000-digit zeros (zeros2). Conventions (derived from
   xi(s) = pi^(-s/2) Gamma(s/2) zeta(s), even test function h):

     Sum_gamma h(gamma) = h(i/2) + h(-i/2)
        + (1/2pi) Int h(r) [Re psi(1/4 + ir/2) - log pi] dr
        - 2 Sum_{p,m} (log p) p^(-m/2) g(m log p),
     g(u) = (1/2pi) Int h(r) e^(-iru) dr.

   Test family: h(r) = H(r-T) + H(r+T), H(x) = sinc(bx)^(2q), band a = 2qb,
   so g(u) = 2 cos(Tu) gHat(u) with gHat(u) = N_{2q}(u/(2b))/(2b) a cardinal
   B-spline vanishing identically for |u| > a = log X  -- the prime side
   cuts off at p^m <= X exactly, no truncation error there.

   Error floors are set by the zero-side tail beyond the 100 known zeros:
   sinc^8 (q=4) predicts ~1e-16, sinc^20 (q=10) predicts ~1e-30; matching
   both floors at three centers T validates every constant in the formula.

   Self-test 0 checks the (H, gHat) Fourier pair numerically first.
*)

prec = 120;
X = 10^5;
a = Log[X];

(* ---- zeros2 parser: blocks separated by blank lines, '.' starts a block *)
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];
lines = ReadList[FileNameJoin[{scriptDir, "..", "refs", "odlyzko", "zeros2"}], String, NullRecords -> True];
zstrs = {}; cur = "";
Do[Module[{l = StringTrim[line]},
   Which[
    l === "", If[cur =!= "", AppendTo[zstrs, cur]; cur = ""],
    StringContainsQ[l, "."], If[cur =!= "", AppendTo[zstrs, cur]]; cur = l,
    True, cur = cur <> l]],
  {line, lines}];
If[cur =!= "", AppendTo[zstrs, cur]];
parseDec[s_String] := Module[{parts = StringSplit[StringTrim[s], "."]},
  FromDigits[parts[[1]]] + FromDigits[parts[[2]]]/10^StringLength[parts[[2]]]];
gammas = parseDec /@ zstrs;
Print[Length[gammas], " zeros parsed; gamma_1 = ", N[gammas[[1]], 30]];
Print["digits per zero (min) = ", Min[StringLength /@ zstrs] - 4];

(* ---- test function family ---- *)
(* cardinal B-spline N_m on [-m/2, m/2]; ~log10 C(m, m/2) digits cancel,
   so evaluate at full working precision *)
bspline[m_, y_] := If[Abs[y] >= m/2, 0,
  Sum[(-1)^k Binomial[m, k] With[{t = y + m/2 - k}, If[t > 0, t^(m - 1), 0]],
    {k, 0, m}]/(m - 1)!];

HH[q_, b_, x_] := If[x == 0, 1, (Sin[b x]/(b x))^(2 q)];        (* also complex x *)
gHat[q_, b_, u_] := bspline[2 q, u/(2 b)]/(2 b);

(* ---- self-test 0: Fourier pair ---- *)
Module[{q = 3, b, u, lhs, rhs},
  b = N[a/6, 40];
  Do[
    lhs = gHat[q, b, N[u, 40]];
    rhs = NIntegrate[HH[q, b, x] Cos[x u], {x, 0, Infinity},
        WorkingPrecision -> 40, PrecisionGoal -> 25, MaxRecursion -> 14]/Pi;
    Print["pair check u=", u, ":  gHat=", N[lhs, 20], "  NInt=", N[rhs, 20],
      "  diff=", N[lhs - rhs, 5]],
    {u, {0, 1, 5, N[a, 40] - 1/2}}]];

(* ---- formula pieces ---- *)
primes = Prime[Range[PrimePi[X]]];

zeroSide[T_, q_, b_] := 2 Sum[
   HH[q, b, N[g - T, prec]] + HH[q, b, N[g + T, prec]], {g, gammas}];

primeSide[T_, q_, b_] := -2 Sum[Module[{lp = Log[N[p, prec]], s = 0, u},
     Do[u = m lp;
       If[u < a, s += lp p^(-m/2) 2 Cos[N[T, prec] u] gHat[q, b, u]],
       {m, 1, Floor[Log[X]/Log[p]]}];
     s], {p, primes}];

poleTerm[T_, q_, b_] := 2 Re[HH[q, b, N[I/2 - T, prec]] + HH[q, b, N[I/2 + T, prec]]];

archTerm[T_, q_, b_, xmax_] := Module[{brks, integrand},
  integrand[x_?NumericQ] :=
    HH[q, b, x] (Re[PolyGamma[0, 1/4 + I (N[T, 80] + x)/2]] - Log[Pi]);
  brks = Table[k Pi/b, {k, 1, Floor[xmax b/Pi]}];
  brks = Sort[Join[-brks, {0}, brks]];
  (1/Pi) NIntegrate[integrand[x],
    Evaluate[Flatten[{x, -xmax, brks, xmax}]],
    WorkingPrecision -> 80, PrecisionGoal -> 40, AccuracyGoal -> 45,
    MaxRecursion -> 10]];

checkIdentity[T_, q_] := Module[{b, xmax, zs, ps, pt, at, rhs, diff},
  b = N[a, prec]/(2 q);
  (* truncate arch integral where H drops below the zero-side tail floor *)
  xmax = Min[174, (10^(34/(2 q)))/N[b]];
  zs = zeroSide[T, q, b];
  ps = primeSide[T, q, b];
  pt = poleTerm[T, q, b];
  at = archTerm[T, q, b, xmax];
  rhs = pt + at + ps;
  diff = zs - rhs;
  Print["T=", N[T, 12], "  q=", q];
  Print["  zeroSide  = ", N[zs, 35]];
  Print["  primeSide = ", N[ps, 35]];
  Print["  archTerm  = ", N[at, 35]];
  Print["  poleTerm  = ", N[pt, 10]];
  Print["  LHS-RHS   = ", N[diff, 6], "   |diff| = 10^", N[Log10[Abs[diff] + 10^-60], 4]];
  diff];

Tlist = {gammas[[50]], 1435/10, 100};
Do[checkIdentity[T, q], {q, {4, 10}}, {T, Tlist}];
