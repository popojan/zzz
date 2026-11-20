#!/usr/bin/env python3
"""How far are the zeros in a `zzz --loop` checkpoint from the true ones?

    ./loop-refzeros 30000 2500 > ref_30000.tsv           # see loop-refzeros.c (FLINT)
    python3 doc/ghy/loop-zero-accuracy.py zzz-loop.state ref_*.tsv

Errors are in units of the local mean gap 2 pi / ln(gamma / 2 pi).  Reference points:
box-smoothed forward (default) ~0.025 rms, plain method B (--loop-boxc 0) ~0.040, a
fabricated zero (pre-2026-11 bracket fallback) ~1.8.  See doc/notes/loop-margins.md.
"""
import math
import sys

import numpy as np

state, refs = sys.argv[1], sys.argv[2:]
zeros, n = [], None
with open(state) as f:
    for line in f:
        if line.startswith("nzeros"):
            n = int(line.split()[1])
            zeros = np.array([float(next(f)) for _ in range(n)])
print(f"{state}: {n} zeros up to {zeros[-1]:.1f}")
print(f"{'index':>9} {'height':>10} | err/gap:  rms    median    p99     max  | worst index")
for ref in sorted(refs, key=lambda r: np.loadtxt(r, max_rows=1)[0]):
    d = np.loadtxt(ref)
    idx, g = d[:, 0].astype(int), d[:, 1]
    ok = idx <= n
    if ok.sum() < 10:
        continue
    idx, g = idx[ok], g[ok]
    e = (zeros[idx - 1] - g) / (2 * math.pi / np.log(g / 2 / math.pi))
    ae = np.abs(e)
    print(f"{idx[0]:9d} {g[0]:10.0f} |         {e.std():6.4f} {np.median(ae):7.4f} {np.quantile(ae, .99):7.4f}"
          f" {ae.max():7.4f} | {idx[ae.argmax()]}")
