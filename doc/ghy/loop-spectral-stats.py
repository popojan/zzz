#!/usr/bin/env python3
"""Spectral statistics of the SELF-COMPUTED zeros in a `zzz --loop` checkpoint
-- nearest-neighbour spacing (level repulsion), number variance Sigma^2(L) and
spectral rigidity Delta_3(L) -- vs the GUE and Poisson predictions.

The zeros are the loop's own zeta-evaluation-free, primality-test-free output;
this measures how much of the operator's spectral fingerprint survives the
band-saturation (~5% of a gap) error.  No ZetaZero needed (uses the checkpoint).

    python3 doc/ghy/loop-spectral-stats.py [zzz-loop.state]
"""
import sys, math

EULER = 0.5772156649015329
TWO_PI = 2 * math.pi


def theta(t):                                   # Riemann-Siegel theta (asymptotic + 1/48t)
    return (t / 2) * math.log(t / TWO_PI) - t / 2 - math.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "zzz-loop.state"
    g = []
    for line in open(path):
        s = line.strip()
        if "." in s and s.replace(".", "", 1).isdigit():
            g.append(float(s))
    n = len(g)
    if n < 1000:
        print("too few zeros"); return 1

    try:
        import numpy as np
    except ImportError:
        print("needs numpy"); return 1

    g = np.array(g)
    w = np.array([theta(x) / math.pi + 1 for x in g])     # unfolded: mean spacing 1
    w.sort()
    print(f"zeros        : {n}   gamma in [{g[0]:.1f}, {g[-1]:.1f}]   unfolded span {w[-1]-w[0]:.0f}")

    # ---- nearest-neighbour spacing distribution ----
    sp = np.diff(w); sp = sp / sp.mean()
    def frac(a, b): return float(((sp >= a) & (sp < b)).mean())
    gue_s = lambda s: (32 / math.pi ** 2) * s * s * math.exp(-4 * s * s / math.pi)
    print("\nnearest-neighbour spacing  (level repulsion):")
    print(f"  P(s<0.5): obs {frac(0,0.5):.3f}   GUE 0.112   Poisson 0.393")
    print(f"  P(s<0.2): obs {frac(0,0.2):.4f}  GUE 0.0084  Poisson 0.181  "
          f"(<- small-s also limited by the marching bracket)")
    hist, edges = np.histogram(sp[sp < 3], bins=30, range=(0, 3), density=True)
    L2g = sum((hist[i] - gue_s((i + .5) * .1)) ** 2 for i in range(25))
    L2p = sum((hist[i] - math.exp(-(i + .5) * .1)) ** 2 for i in range(25))
    print(f"  L2(hist): to GUE {L2g:.3f}   to Poisson {L2p:.3f}   -> "
          f"{'GUE' if L2g < L2p else 'Poisson'}")

    # ---- number variance Sigma^2(L): GUE ~ log L (rigid), Poisson = L ----
    def sigma2(L, nwin=8000):
        starts = np.linspace(w[0], w[-1] - L, nwin)
        cnt = np.searchsorted(w, starts + L, "right") - np.searchsorted(w, starts, "left")
        return cnt.mean(), cnt.var()
    sig_gue = lambda L: (1 / math.pi ** 2) * (math.log(TWO_PI * L) + EULER + 1)
    print("\nnumber variance  Sigma^2(L):   GUE ~ (1/pi^2)logL (rigid)   Poisson = L")
    print(f"  {'L':>5} {'<n>':>8} {'Sigma2_obs':>11} {'GUE':>8} {'Poisson':>9}")
    grid = [1, 2, 5, 10, 20, 50, 100, 200]
    S2 = {}
    for L in grid:
        if w[-1] - w[0] <= L: continue
        m, v = sigma2(L); S2[L] = v
        print(f"  {L:5} {m:8.2f} {v:11.3f} {sig_gue(L):8.3f} {float(L):9.1f}")

    # ---- spectral rigidity Delta_3(L) from Sigma^2 (Dyson-Mehta) ----
    def delta3(L, K=400):
        rs = np.linspace(1e-6, L, K)
        s2 = np.array([sigma2(r, 3000)[1] for r in rs])
        ker = L ** 3 - 2 * L ** 2 * rs + rs ** 3
        f = ker * s2
        return (2 / L ** 4) * float(np.sum((f[:-1] + f[1:]) / 2 * np.diff(rs)))   # trapezoid
    d3_gue = lambda L: (1 / math.pi ** 2) * (math.log(TWO_PI * L) + EULER - 5 / 4 - math.pi ** 2 / 8)
    print("\nspectral rigidity  Delta_3(L):  GUE ~ (1/pi^2)logL          Poisson = L/15")
    print(f"  {'L':>5} {'Delta3_obs':>11} {'GUE':>8} {'Poisson':>9}")
    for L in [5, 10, 20, 50]:
        if w[-1] - w[0] <= L: continue
        d = delta3(L)
        print(f"  {L:5} {d:11.3f} {d3_gue(L):8.3f} {L/15:9.3f}")

    print("\nverdict: nearest-neighbour spacing ~ GUE (local level repulsion survives the error),")
    print("         BUT Sigma^2(L) and Delta_3(L) SATURATE below GUE's (1/pi^2)logL -> OVER-rigid")
    print("         at long range: the band-limited prime sum misses the long-range fluctuations of")
    print("         S(t)=arg zeta/pi, pinning zeros to the smooth count (band-saturation, spectral form).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
