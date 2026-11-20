#!/usr/bin/env python3
"""Decision margins of `zzz --loop`, from its own decision log.

    ./zzz --loop --loop-fresh --loop-declog run.dec > run.primes     # Ctrl+C when bored
    python3 doc/ghy/loop-margins.py run.dec [bins]

Each log line is one backward decision:
    iter  x  nz  gamma_top  gate_frac  median  code  known_power  r_0 .. r_5
code: P prime-height, C clear composite, a ambiguous but already known, A ambiguous
(scan stops), G out of reach (scan stops).  The truth labels come from a sieve here;
the loop never sees them.

Per x-range the table gives the commit lag a = gamma_top / (2 pi x), the median statistic
for primes and for composites, the distance of the composite cloud from the 0.40 commit
threshold in its own standard deviations, and var * a * ln x, which is flat if the noise
shrinks like 1/sqrt(a ln x).  alpha = z^2 / (2 ln x) is the exponent of the per-candidate
false-positive probability x^-alpha under a Gaussian tail; alpha > 1 is summable.
See doc/notes/loop-margins.md.
"""
import math
import sys

import numpy as np

TAU_C, TAU_LO = 0.40, 0.25


def main():
    path = sys.argv[1]
    nbins = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    dec, stops = [], []
    for line in open(path):
        p = line.split("\t")
        if len(p) < 14:
            continue
        row = (int(p[0]), int(p[1]), float(p[3]), float(p[5]), p[6], int(p[7]))
        (stops if p[6] in "GA" else dec).append(row)
    x = np.array([r[1] for r in dec])
    _, first = np.unique(x, return_index=True)          # one decision per x
    dec = [dec[i] for i in first]
    x = np.array([r[1] for r in dec])
    top = np.array([r[2] for r in dec])
    m = np.array([r[3] for r in dec])
    code = np.array([r[4] for r in dec])
    kp = np.array([r[5] for r in dec])

    n = int(x.max()) + 1
    sv = np.ones(n + 1, bool)
    sv[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    isp = sv[x]
    a = top / (2 * math.pi * x)

    print(f"{len(x)} decisions, x = {x.min()} .. {x.max()}")
    print("committed as prime but composite:", int(((code == 'P') & ~isp & (kp == 0)).sum()),
          "  declared composite but prime:", int(((code == 'C') & isp).sum()))
    lo = max(int(x.min()), 600)                          # below ~600 the lag is still settling
    edges = np.unique(np.round(np.exp(np.linspace(math.log(lo), math.log(x.max() + 1), nbins + 1))).astype(int))
    print(f"{'x range':>15} {'ln x':>5} {'lag a':>5} | primes    n   mean     sd    min |"
          f" composites  n    mean     sd    max | to {TAU_C:.2f}  var*a*lnx  alpha")
    zs = []
    for e0, e1 in zip(edges[:-1], edges[1:]):
        s = (x >= e0) & (x < e1)
        pr, co = s & isp, s & ~isp & (kp == 0)           # prime powers of known primes are bookkeeping
        if pr.sum() < 5 or co.sum() < 20:
            continue
        lx, am = math.log((e0 + e1) / 2), float(np.median(a[s]))
        mp, mc = m[pr], m[co]
        z = (TAU_C - mc.mean()) / mc.std()
        zs.append((mc - mc.mean()) / mc.std())
        print(f"{e0:7d}-{e1:7d} {lx:5.2f} {am:5.2f} | {len(mp):10d} {mp.mean():6.3f} {mp.std():6.3f} {mp.min():6.3f} |"
              f" {len(mc):9d} {mc.mean():+7.3f} {mc.std():6.4f} {mc.max():6.3f} | {z:5.1f}s"
              f"   {mc.var() * am * lx:7.4f}  {z * z / 2 / lx:5.2f}")
    zs = np.concatenate(zs)
    tail = "  ".join(f"P(z>{t})={(zs > t).mean():.1e} [gauss {0.5 * math.erfc(t / 2 ** .5):.1e}]" for t in (2, 2.5, 3, 3.5))
    print(f"composite tail, standardised per bin, n={len(zs)}: {tail}  max z={zs.max():.2f}")
    q = np.quantile(a[x >= lo], [0, 0.01, 0.5, 0.99, 1])
    print("commit lag a: min %.3f  1%% %.3f  median %.3f  99%% %.3f  max %.3f   (floor ~0.55)" % tuple(q))
    kinds = [r[4] for r in stops]
    print(f"scan stops: {kinds.count('G')} out of reach, {kinds.count('A')} ambiguous")


if __name__ == "__main__":
    sys.exit(main())
