\\ check-loop.gp -- verify the primes in a `zzz --loop` checkpoint are CORRECT
\\ (all actually prime) and COMPLETE (no prime up to the largest one is missing).
\\ This is the external, post-hoc check the loop itself deliberately never does.
\\
\\ Usage:   gp -q doc/ghy/check-loop.gp                 # checks ./zzz-loop.state
\\          STATE=/path/to.state gp -q doc/ghy/check-loop.gp
\\ In the file: prime entries are the bare-integer lines (header lines carry
\\ letters, zero ordinates carry a '.'), so they are trivial to pick out.

checkloop(fname) =
{
  my(L = readstr(fname), P = List(), xk = -1);
  for(i = 1, #L,
    my(v = Vecsmall(L[i]), dig = (#v > 0));
    for(j = 1, #v, if(v[j] < 48 || v[j] > 57, dig = 0; break));
    if(dig, listput(P, eval(L[i])));                       \\ bare integer -> a prime entry
    if(#v >= 6 && v[1] == 120 && v[2] == 107,              \\ "xknown ..." line
       xk = eval(Strchr(Vecsmall(select(c -> c >= 48 && c <= 57, Vec(v))))))
  );
  P = vecsort(Vec(P));
  my(np = #P);
  if(np == 0, print("no primes in ", fname); return);
  my(maxP = P[np], tru = primes([2, maxP]));
  my(wrong   = select(p -> !isprime(p), P),               \\ listed but composite
     missing = setminus(Set(tru), Set(P)),                \\ prime <= maxP not listed
     dups    = np - #Set(P));
  print("file        : ", fname);
  print("listed primes: ", np, "   (largest ", maxP, ")");
  print("X_known     : ", if(xk >= 0, xk, "?"));
  print("not prime   : ", if(#wrong,   wrong,   "none"));
  print("missing     : ", if(#missing, missing, "none"));
  print("duplicates  : ", dups);
  if(xk >= 0 && primepi(xk) != np,
     print("note        : pi(X_known)=", primepi(xk), " != ", np,
           " (primes between ", maxP, " and X_known=", xk, " not all found)"));
  if(#wrong == 0 && #missing == 0 && dups == 0,
     print("RESULT      : CORRECT and COMPLETE -- exactly the primes up to ", maxP),
     print("RESULT      : *** PROBLEM ***"));
}

{
  my(f = getenv("STATE"));
  if(f == "", f = "zzz-loop.state");
  checkloop(f);
}
