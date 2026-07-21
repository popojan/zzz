#!/usr/bin/env python3
"""Verify the primes in a `zzz --loop` checkpoint are CORRECT (every listed
number is prime) and COMPLETE (no prime up to X_known is missing) -- the
external, post-hoc check the loop itself deliberately never does.

Pure Python 3 stdlib (no sympy), streams the file (handles huge checkpoints).
Prime entries are the bare-integer lines; header lines carry letters, zero
ordinates carry a '.'.

    python3 doc/ghy/check-loop.py [zzz-loop.state]
"""
import sys


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "zzz-loop.state"
    primes, xknown = [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.isdigit():              # bare integer (no '.') -> a prime entry
                primes.append(int(s))
            elif s.startswith("xknown"):
                xknown = int(s.split()[1])
            # other header lines (letters) and zeros (contain '.') are skipped
    if not primes:
        print(f"no primes in {path}")
        return 1

    pset = set(primes)
    mx = max(primes)
    bound = max(mx, xknown or 0)

    sieve = bytearray([1]) * (bound + 1)         # 1 = prime
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(bound ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, bound + 1, i)))

    wrong = sorted(p for p in pset if not sieve[p])               # listed but composite
    missing = [p for p in range(2, bound + 1) if sieve[p] and p not in pset]
    dups = len(primes) - len(pset)

    def show(xs):
        return "none" if not xs else (str(xs[:10]) + (" ..." if len(xs) > 10 else ""))

    print(f"file         : {path}")
    print(f"listed primes: {len(primes)}  (largest {mx})")
    print(f"X_known      : {xknown if xknown is not None else '?'}")
    print(f"not prime    : {show(wrong)}")
    print(f"missing      : {show(missing)}   (primes <= {bound} not listed)")
    print(f"duplicates   : {dups}")
    ok = not wrong and not missing and not dups
    print("RESULT       : " +
          (f"CORRECT and COMPLETE -- exactly the primes up to {bound}"
           if ok else "*** PROBLEM ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
