#!/usr/bin/env python3
"""What do zeros located from the primes <= X know about the integers beyond X?

    ./loop-fixed-cutoff 100 120 200000 0.375 z100.bin        # see loop-fixed-cutoff.c
    python3 doc/ghy/loop-fixed-cutoff.py 100 z100.bin

The zeros are binned by the forward margin they were located at,
kappa = 2 pi ln X / ln(gamma / 2 pi), and for several classes of integers x the table shows
the coherence of each bin: -sum_k cos(gamma_k ln x), in units of the peak the TRUE zeros of
that bin would give a prime at x (Landau: (d gamma / 2 pi) ln x / sqrt x).  So 1 = "reads as a
prime", 0 = "reads as a composite", and the true zeros would give 1/m at a prime power p^m.

Result (X = 100): primes beyond X read 0.5-1.0 at every margin.  Composites with a prime
factor <= X read 0.  Composites p*q with BOTH factors > X read ~0 while kappa >~ 5, but
0.5-0.7 once kappa <~ 4, and p^2 with p > X reads like a prime rather than 1/2: near
kappa = pi the forward step only transmits "x has no prime factor <= X" -- the sieve of
Eratosthenes, and X = sqrt(T/2pi) is where that stops identifying primes (x < X^2).
See doc/notes/loop-margins.md.
"""
import math
import random
import sys

import numpy as np

X = int(sys.argv[1])
Z = np.fromfile(sys.argv[2])
KAPPA = [9, 2 * math.pi, 5.5, 5, 4.5, 4, 3.5, math.pi, 2.9, 2.7]
edges = [2 * math.pi * X ** (2 * math.pi / k) for k in KAPPA]
idx = np.searchsorted(Z, edges)

TOP = 5 * X * X
spf = np.zeros(TOP + 1, np.int64)                    # smallest prime factor
for i in range(2, TOP + 1):
    if spf[i] == 0:
        spf[i::i] = np.where(spf[i::i] == 0, i, spf[i::i])


def is_prime(n): return spf[n] == n
def is_semiprime(n): return not is_prime(n) and is_prime(n // spf[n])


def profile(xs):
    acc = np.zeros(len(KAPPA) - 1)
    for x in xs:
        lx = math.log(x)
        cs = np.concatenate([[0], np.cumsum(np.cos(Z * lx))])
        for b in range(len(KAPPA) - 1):
            if idx[b + 1] > idx[b]:
                dT = Z[min(idx[b + 1], len(Z)) - 1] - Z[idx[b]]
                acc[b] += -(cs[idx[b + 1]] - cs[idx[b]]) / (dT / (2 * math.pi) * lx / math.sqrt(x))
    return acc / len(xs)


def show(label, xs, cap=150):
    xs = list(xs)[:cap]
    if xs:
        print(f"{label:<40} n={len(xs):4d} : " + " ".join(f"{v:6.2f}" for v in profile(xs)))


primes = [n for n in range(2, TOP + 1) if is_prime(n)]
print(f"X={X}: {len(Z)} zeros, gamma = {Z[0]:.0f} .. {Z[-1]:.0f}")
print("kappa bins: " + "  ".join(f"{a:.2f}-{b:.2f}" for a, b in zip(KAPPA[1:], KAPPA[:-1])) + "   (left = highest margin)")
print(f"{'':<51}" + " ".join(f"{k:6.2f}" for k in KAPPA[1:]))
show("primes <= X (in the forward sum)", [p for p in primes if p <= X])
for lo, hi in [(X, 2 * X), (2 * X, 5 * X), (5 * X, 20 * X), (20 * X, X * X), (X * X, 2 * X * X), (2 * X * X, 5 * X * X)]:
    show(f"primes in ({lo},{hi}]", [p for p in primes if lo < p <= hi])
show("p^2 > X, p <= X       (true zeros: 0.5)", [p * p for p in primes if p <= X < p * p])
show("p^2, p > X            (true zeros: 0.5)", [p * p for p in primes if X < p and p * p <= TOP])
show("p*q, both > X         (true zeros: 0)", [n for n in range(X * X, TOP) if spf[n] > X and is_semiprime(n)])
random.seed(1)
show("composites > X^2 with a factor <= X", [n for n in random.sample(range(X * X, TOP), 600) if spf[n] <= X and not is_prime(n)])
show("composites in (2X, X^2)", [n for n in random.sample(range(2 * X, X * X), 600) if not is_prime(n)])
