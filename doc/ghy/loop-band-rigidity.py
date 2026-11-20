#!/usr/bin/env python3
"""Band-edge -> spectral-crossover scaling law.

Recompute one fixed block of zeros (#A..#B, fixed height) by method B at a SWEEP
of prime bands X, and measure the number variance Sigma^2(L) of each, vs the TRUE
zeros for the same block (X -> infinity).  Tests whether the band sets the
long-range rigidity:  prediction Sigma^2_sat(X) ~ (1/pi^2) log log X, climbing to
the true value as X -> Riemann-Siegel scale sqrt(T/2pi).

Needs the true block once (Mathematica):
    wolframscript -code 'Export["/tmp/true_block.txt",
        Table[N[Im[ZetaZero[k]]], {k, 30000, 31800}], "Table"]'
    python3 doc/ghy/loop-band-rigidity.py            # writes loop-band-rigidity.png
See doc/notes/spectral-shadow.md.
"""
import math, sys
import numpy as np

PI, TWO_PI = math.pi, 2 * math.pi
A, B = 30000, 31800
Xs = [32, 48, 64, 96, 128, 256, 512, 1024]
Ls = [2, 3, 5, 8, 12, 18, 25]

def N0(t):  u = t / TWO_PI; return u * math.log(u / math.e) + 0.875
def dN0(t): return math.log(t / TWO_PI) / TWO_PI
def theta(t): return (t / 2) * math.log(t / TWO_PI) - t / 2 - PI / 8 + 1 / (48 * t)
def unfold(g): return np.sort(np.array([theta(x) / PI + 1 for x in g]))

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def prime_powers(X):                                   # (m log p, (1/m) p^{-m/2}) for p^m<=X
    fr, am = [], []
    for p in sieve(X):
        m, pm = 1, int(p)
        while pm <= X:
            fr.append(m * math.log(p)); am.append((1.0 / m) * p ** (-0.5 * m)); m += 1; pm *= p
    return np.array(fr), np.array(am)

def locate_block(X):                                   # method-B zeros #A..#B at band X
    fr, am = prime_powers(X)
    FB = lambda t: N0(t) - float(am @ np.sin(t * fr)) / PI
    t = 25000.0
    for _ in range(80): t -= (N0(t) - (A - 0.5)) / dN0(t)
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

def sigma2(w, L, nwin=4000):
    s = np.linspace(w[0], w[-1] - L, nwin)
    return float((np.searchsorted(w, s + L, "right") - np.searchsorted(w, s, "left")).var())

true = [float(x) for x in open("/tmp/true_block.txt") if x.strip()][:B - A + 1]
T = sum(true) / len(true)
Xrs = math.sqrt(T / TWO_PI)
wT = unfold(true)
sig_gue = lambda L: (1 / PI ** 2) * (math.log(TWO_PI * L) + 0.5772156649 + 1)
llX = lambda X: (1 / PI ** 2) * math.log(math.log(X))     # ~ saturation prediction scale

gue_s = lambda s: (32 / PI ** 2) * s * s * math.exp(-4 * s * s / PI)
def shortrange(w):                                     # local spacing: P(s<0.5), L2-to-GUE
    sp = np.diff(w); sp = sp / sp.mean()
    h, _ = np.histogram(sp[sp < 3], bins=30, range=(0, 3), density=True)
    return float((sp < 0.5).mean()), sum((h[i] - gue_s((i + .5) * .1)) ** 2 for i in range(25))

curves, short = {}, {}
for X in Xs:
    w = unfold(locate_block(X))
    curves[X] = [sigma2(w, L) for L in Ls]; short[X] = shortrange(w)
trueC = [sigma2(wT, L) for L in Ls]; trueSR = shortrange(wT)

print(f"block #{A}-{B}: gamma~{T:.0f}, Riemann-Siegel band X_RS=sqrt(T/2pi)~{Xrs:.0f}\n")
hdr = "  L  " + "".join(f"X={X:<5}" for X in Xs) + "  TRUE   GUE"
print(hdr)
for i, L in enumerate(Ls):
    print(f"{L:4} " + "".join(f"{curves[X][i]:6.3f} " for X in Xs)
          + f"{trueC[i]:6.3f} {sig_gue(L):6.3f}")
print("\nsaturation value (Sigma^2 at L=%d) vs band:" % Ls[-1])
print("   X     Sig2_sat   (1/pi^2)loglogX")
for X in Xs:
    print(f"  {X:5} {curves[X][-1]:8.3f}   {llX(X):10.3f}")
print(f"  TRUE  {trueC[-1]:8.3f}   {(1/PI**2)*math.log(math.log(T)):10.3f}  (loglog T)")

print(f"\nSHORT range (local spacing) vs band:   kappa = gap*logX, threshold pi at X_RS~{Xrs:.0f}")
print("   X    kappa   P(s<0.5)   L2-to-GUE   (GUE: 0.112 / 0.0 ;  Poisson P=0.393)")
gap = TWO_PI / math.log(T / TWO_PI)
for X in Xs:
    print(f"  {X:5} {gap*math.log(X):6.2f}   {short[X][0]:7.3f}    {short[X][1]:8.3f}")
print(f"  TRUE         {trueSR[0]:7.3f}    {trueSR[1]:8.3f}")

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for X in Xs: ax[0].plot(Ls, curves[X], "-o", ms=3, label=f"X={X}")
    ax[0].plot(Ls, trueC, "k--", lw=2, label="TRUE")
    ax[0].plot(Ls, [sig_gue(L) for L in Ls], "0.6", ls=":", label="GUE logL")
    ax[0].axvline(0); ax[0].set_xlabel("L"); ax[0].set_ylabel(r"$\Sigma^2(L)$")
    ax[0].set_title(f"number variance vs band  (gamma~{T:.0f}, X_RS~{Xrs:.0f})"); ax[0].legend(fontsize=7)
    ax[1].semilogx(Xs, [short[X][0] for X in Xs], "s-", label="P(s<0.5)")
    ax[1].axhline(0.112, color="b", ls=":", label="GUE 0.112")
    ax[1].axhline(0.393, color="r", ls=":", label="Poisson 0.393")
    ax[1].axhline(trueSR[0], color="k", ls="--", label="TRUE")
    ax[1].axvline(Xrs, color="g", ls=":", label=r"$X_{RS}$ ($\kappa=\pi$)")
    ax[1].set_xlabel("band X"); ax[1].set_ylabel("P(s<0.5)  (level repulsion)"); ax[1].legend(fontsize=7)
    ax[1].set_title("local spacing onset vs band edge")
    fig.tight_layout(); fig.savefig("loop-band-rigidity.png", dpi=110)
    print("\nwrote loop-band-rigidity.png")
except ImportError:
    print("\n(matplotlib absent; table only)")
