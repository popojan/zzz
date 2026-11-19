#!/usr/bin/env python3
"""Replay the `zzz --loop` backward decision rule at heights the loop cannot reach.

The statistic the loop thresholds is modelled, per truncation b, as

    r_b(x) = A_b * ( -1 + 2 (u_b/pi) C  sum_n Lambda(n) sinc(u_b (x - n)) ) + jitter,

u_b = gamma(FR_b * NZ) / x, A_b = pi / (2 u_b ln x): the explicit-formula peak of each
prime, smeared by the hard truncation of the zero sum, at coherence C (how strongly the
loop's zeros carry a prime that was not in the forward sum yet; ~0.90 box-B, ~0.72
plain-B), plus prefix-sum noise of variance c / (a_b ln x).  The primes are the real ones
around x0 (segmented sieve), so twin primes and denser clusters appear with their true
frequency.  Against the loop's own decision log the model reproduces the per-truncation
composite readings with correlation ~0.8 and the residual gives c ~ 0.021.  (In these units
a prime reads C - 1/(4 a_b ln x) and an isolated composite -1/(4 a_b ln x).)

For each lag a = gamma_top / (2 pi x) the table shows how often the statistical reach gate
passes, and what the median rule would then decide.  The point: false primes (composites
between twin primes) appear for a <~ 0.55 at EVERY height, none for a >= 0.6; the loop
commits at a ~ 0.9 and the explicit floor U_MIN = 2.5 forbids a < ~0.75.

    python3 doc/ghy/loop-replay.py X0 [SEGMENT=1e6] [C=0.90] [jitter c=0.021]
    for x in 6e3 2.5e5 1e7 1e9 1e12 1e15; do python3 doc/ghy/loop-replay.py $x; done

See doc/notes/loop-margins.md.  Needs numpy + scipy.
"""
import math
import sys

import numpy as np
from scipy.signal import fftconvolve

FR = [.5, .6, .7, .8, .9, 1.0]
LAGS = [1.3, 1.1, 1.0, 0.9, 0.8, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4]
TAU_C, TAU_LO, TAULO_GATE, FRAC_OK, GW, U_MIN = 0.40, 0.25, 0.12, 0.70, 30, 2.5
W = 3000                                             # kernel half-width (candidates)


def N0(t):
    u = t / (2 * math.pi)
    return u * math.log(u / math.e) + 0.875


def n0inv(n):
    t = 2 * math.pi * n / math.log(n + 3) + 10
    for _ in range(80):
        t -= (N0(t) - n) / (math.log(t / 2 / math.pi) / (2 * math.pi))
    return t


def base_primes(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for a in range(2, int(n ** .5) + 1):
        if s[a]:
            s[a * a::a] = False
    return np.nonzero(s)[0].astype(np.int64)


def segment(x0, length, bp):
    """primality of x0 .. x0+length-1 (x0 well above sqrt of the top)"""
    seg = np.ones(length, bool)
    for p in bp[bp <= length]:
        seg[(-x0) % int(p)::p] = False
    big = bp[bp > length]
    if len(big):
        st = (-np.int64(x0)) % big
        seg[st[st < length]] = False
    for p in bp[(bp >= x0) & (bp < x0 + length)]:    # small x0: a base prime crossed itself out
        seg[p - x0] = True
    return seg


def main():
    x0 = int(float(sys.argv[1]))
    L = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 6
    C = float(sys.argv[3]) if len(sys.argv) > 3 else 0.90
    cj = float(sys.argv[4]) if len(sys.argv) > 4 else 0.021
    L = min(L, x0)
    seg = segment(x0 - W, L + 2 * W, base_primes(int((x0 + L + W) ** .5) + 2))
    lx = math.log(x0)
    lam = seg.astype(float) * lx
    isp = seg[W:W + L]
    comp = ~isp
    twin = np.zeros(L, bool)
    twin[1:-1] = comp[1:-1] & isp[:-2] & isp[2:]
    rng = np.random.default_rng(7)
    dd = np.arange(-W, W + 1)
    print(f"x0={x0:.3g}  ln x={lx:.2f}  {L} candidates, {isp.sum()} primes, {twin.sum()} twin-sandwiched composites"
          f"  C={C}  c={cj}")
    print(" lag a  g(.5NZ)/x  gate pass% | of those: FALSE PRIME  missed  ambiguous | all x: false  missed |"
          " composite median max   prime median min | twin-sandwich median mean / max")
    for a in LAGS:
        T = 2 * math.pi * a * x0
        N = N0(T)
        Tb = [n0inv(fr * N) if fr < 1 else T for fr in FR]
        R = np.empty((6, L))
        walk, prevT = np.zeros(L), 0.0
        for b in range(6):
            u = Tb[b] / x0
            A = math.pi / (2 * u * lx)                       # loop.c's scale (2 gamma/pi) ln x / x
            conv = fftconvolve(lam, np.sinc(u * dd / math.pi), mode="valid")
            walk = walk + math.sqrt(Tb[b] - prevT) * rng.standard_normal(L)      # prefix sums share noise
            prevT = Tb[b]
            ab = Tb[b] / (2 * math.pi * x0)
            R[b] = A * (-1 + 2 * (u / math.pi) * C * conv) + walk / math.sqrt(Tb[b]) * math.sqrt(cj / (ab * lx))
        bad = (R[0] >= TAULO_GATE).astype(float)
        cs = np.concatenate([[0], np.cumsum(bad)])
        lo, hi = np.clip(np.arange(L) - GW, 0, L), np.clip(np.arange(L) + GW + 1, 0, L)
        gate = 1 - (cs[hi] - cs[lo]) / (hi - lo) >= FRAC_OK
        srt = np.sort(R, axis=0)
        med = 0.5 * (srt[2] + srt[3])
        fp, miss = comp & (med > TAU_C), isp & (med < TAU_LO)
        amb = (med >= TAU_LO) & (med <= TAU_C)
        flag = "" if Tb[0] / x0 >= U_MIN else "   <- forbidden by U_MIN"
        print(f"{a:6.2f} {Tb[0] / x0:9.2f} {100 * gate.mean():10.3f}  | {(gate & fp).sum():20d} {(gate & miss).sum():7d}"
              f" {(gate & amb).sum():10d} | {fp.sum():12d} {miss.sum():7d} | {med[comp].max():20.3f}"
              f" {med[isp].min():18.3f} | {med[twin].mean():14.3f} / {med[twin].max():.3f}{flag}")


if __name__ == "__main__":
    sys.exit(main())
