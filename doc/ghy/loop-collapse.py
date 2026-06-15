#!/usr/bin/env python3
"""Universality of the band-saturation transition.

At several heights T, recompute a block of zeros by method B over a sweep of bands
X, and plot the spectral order parameter P(s<0.5) (level repulsion) against
kappa = gap*logX.  If the curves from all heights COLLAPSE onto one curve, the
crystal->GUE transition is universal in kappa, located at kappa=pi -- i.e.
X=sqrt(T/2pi) is a scaling LAW, and the information / spectral / cost faces of the
wall are one universal crossover.  Pure method-B (no ZetaZero).  -> loop-collapse.png
See doc/notes/spectral-shadow.md.
"""
import math
import numpy as np

PI, TWO_PI = math.pi, 2 * math.pi
HEIGHTS = [5000, 25000, 60000, 120000]          # block start indices
BLOCK = 2500
KAPPAS = [2.0, 2.5, 2.8, 3.14, 3.5, 4.0, 4.5, 5.0, 6.0]

def N0(t):  u = t / TWO_PI; return u * math.log(u / math.e) + 0.875
def dN0(t): return math.log(t / TWO_PI) / TWO_PI
def theta(t): return (t / 2) * math.log(t / TWO_PI) - t / 2 - PI / 8 + 1 / (48 * t)
def unfold(g): return np.sort(np.array([theta(x) / PI + 1 for x in g]))

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def prime_powers(X):
    fr, am = [], []
    for p in sieve(int(X)):
        m, pm = 1, int(p)
        while pm <= X:
            fr.append(m * math.log(p)); am.append((1.0 / m) * p ** (-0.5 * m)); m += 1; pm *= p
    return np.array(fr), np.array(am)

def locate_block(A, B, X):
    fr, am = prime_powers(X)
    FB = lambda t: N0(t) - float(am @ np.sin(t * fr)) / PI
    t = max(20.0, TWO_PI * (A / math.log(A)) ** 0)  # rough; refine by Newton
    t = float(A)                                    # seed; Newton converges fast
    for _ in range(100): t -= (N0(t) - (A - 0.5)) / dN0(t)
    z, prev = [], t
    for k in range(A, B + 1):
        g = TWO_PI / math.log(max(prev, 10) / TWO_PI); tgt = k - 0.5
        lo, hi = prev + 0.25 * g, prev + 1.9 * g
        flo, fhi = FB(lo) - tgt, FB(hi) - tgt
        if flo * fhi > 0:
            lo, hi = prev + 0.02 * g, prev + 2.8 * g; flo, fhi = FB(lo) - tgt, FB(hi) - tgt
        if flo * fhi > 0:
            zk = prev + g
        else:
            for _ in range(48):
                mid = .5 * (lo + hi); fm = FB(mid) - tgt
                if flo * fm <= 0: hi = mid
                else: lo, flo = mid, fm
            zk = .5 * (lo + hi)
        z.append(zk); prev = zk
    return z

def Pls(zeros):
    w = unfold(zeros); sp = np.diff(w); sp = sp / sp.mean()
    return float((sp < 0.5).mean())

curves = {}
print(f"  block of {BLOCK} zeros at each height; sweep band so kappa=gap*logX matches\n")
for A in HEIGHTS:
    t = float(A)
    for _ in range(100): t -= (N0(t) - (A - 0.5)) / dN0(t)
    gap = TWO_PI / math.log(t / TWO_PI)
    row = []
    for kap in KAPPAS:
        X = max(2, int(round(math.exp(kap / gap))))
        row.append((kap, X, Pls(locate_block(A, A + BLOCK, X))))
    curves[A] = (t, gap, row)
    print(f"gamma~{t:7.0f}  gap={gap:.3f}   " +
          "  ".join(f"k{kap:.2f}:P={p:.3f}(X={X})" for kap, X, p in row))

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    for A, (t, gap, row) in curves.items():
        ks = [r[0] for r in row]; ps = [r[2] for r in row]
        ax.plot(ks, ps, "-o", ms=4, label=f"$\\gamma\\sim${t:.0f}")
    ax.axvline(PI, color="g", ls="--", label=r"$\kappa=\pi$")
    ax.axhline(0.112, color="b", ls=":", label="GUE 0.112")
    ax.axhline(0.393, color="r", ls=":", label="Poisson 0.393")
    ax.set_xlabel(r"$\kappa = \mathrm{gap}\cdot\log X$")
    ax.set_ylabel("P(s<0.5)  (level repulsion)")
    ax.set_title("band-saturation transition: collapse in $\\kappa$ across heights")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("loop-collapse.png", dpi=120)
    print("\nwrote loop-collapse.png")
except ImportError:
    print("\n(matplotlib absent; table only)")
