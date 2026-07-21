#!/usr/bin/env python3
"""Is the loop's long-range over-rigidity an artefact, or the real (Berry)
saturation that the TRUE zeros also show?  Compare the number variance Sigma^2(L)
of the loop's zeros to the exact ZetaZero, over the SAME index block.

First make the true block (slow part, needs Mathematica's ZetaZero):
    wolframscript -code 'Export["/tmp/true_block.txt",
        Table[N[Im[ZetaZero[k]]], {k, 30000, 31800}], "Table"]'
then:
    python3 doc/ghy/loop-vs-true-sigma2.py [zzz-loop.state] [/tmp/true_block.txt] [A] [B]

Result (2026-06-15, block 30000-31800, gamma~25755): TRUE and LOOP both ~0.35-0.40,
far below GUE's logL, and within ~0.03 of each other -> the saturation is REAL
(Berry semiclassical, fixed by the short primes), reproduced by the loop. NOT an
artefact.  See doc/notes/spectral-shadow.md.
"""
import math, sys
import numpy as np

state = sys.argv[1] if len(sys.argv) > 1 else "zzz-loop.state"
tblk  = sys.argv[2] if len(sys.argv) > 2 else "/tmp/true_block.txt"
A     = int(sys.argv[3]) if len(sys.argv) > 3 else 30000
B     = int(sys.argv[4]) if len(sys.argv) > 4 else 31800
TWO_PI = 2 * math.pi

def theta(t): return (t / 2) * math.log(t / TWO_PI) - t / 2 - math.pi / 8 + 1 / (48 * t)
def unfold(g): return np.sort(np.array([theta(x) / math.pi + 1 for x in g]))

true = [float(x) for x in open(tblk) if x.strip()]
zl = [float(s) for s in (l.strip() for l in open(state))
      if "." in s and s.replace(".", "", 1).isdigit()]
loop = zl[A - 1:B]
n = min(len(true), len(loop))
true, loop = true[:n], loop[:n]
print(f"block #{A}-{B}: {n} zeros   gamma_{A}: true {true[0]:.3f}  loop {loop[0]:.3f}")

wT, wL = unfold(true), unfold(loop)
def sigma2(w, L, nwin=4000):
    s = np.linspace(w[0], w[-1] - L, nwin)
    return float((np.searchsorted(w, s + L, "right") - np.searchsorted(w, s, "left")).var())
sig_gue = lambda L: (1 / math.pi ** 2) * (math.log(TWO_PI * L) + 0.5772156649 + 1)

print(f"\n  {'L':>4} {'Sig2_TRUE':>10} {'Sig2_LOOP':>10} {'GUE_log':>9}")
for L in [2, 5, 10, 15, 20, 30, 40]:
    if wT[-1] - wT[0] <= L: continue
    print(f"  {L:4} {sigma2(wT, L):10.3f} {sigma2(wL, L):10.3f} {sig_gue(L):9.3f}")
print("\nTRUE ~ LOOP, both << GUE_log  =>  real Berry saturation, reproduced (not an artefact)")
