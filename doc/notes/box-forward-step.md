# Box-averaged forward step: free sinc weights sharpen the loop's method B

*Synthesis of the 2026-07-09/10 experiment. Origin: the orbit session
`2026-07-09-zeta-count-smoothed-ceilings` §9 found that box-averaging (moving
average, window `w`) the truncated critical-line log sum gives a dilogarithm
closed form and a real error reduction against `log ζ(1/2+it)`. This note maps
that result onto zzz, bounds where it can and cannot help, and records the
measured in-loop gain. Gate-check script: `doc/ghy/box-forward-gatecheck.wls`;
implementation: `loop.c::FB` + `--loop-boxc` (default 0.375).*

## The identity

Box-averaging method B's counting function over `[t-w, t+w]` needs no dilog for
the truncated Dirichlet sum — it is exactly a per-term reweighting:

$$\frac{1}{2w}\int_{t-w}^{t+w} F_B(x)\,dx \;=\; N_0 \;+\; \frac1\pi\sum_{p^m\le X}
\frac{p^{-m/2}}{m}\,\operatorname{sinc}(m\,w\log p)\,\sin(m\,t\log p) \;+\; O\!\Big(\frac{w^2}{t}\Big)$$

one multiplicative weight per prime power. Two structural facts:

- **The bisection target is unchanged.** The box average of the *true* staircase
  `N(t)` still crosses `k − 1/2` exactly at `γ_k` whenever the window contains a
  single zero — the smoothing is unbiased for crossing location; bias enters
  only through close pairs when `w` is too large.
- **It is a band-edge taper of the forward prime sum**: `sinc(w log p) ≈ 1` for
  small primes, damping concentrated on the noisy truncation edge — exactly the
  Si/Gibbs neighbour tails that the cutoff-sweep neighbour model identifies as
  B's dominant above-threshold error (corr 0.96–0.98).

The window scales with the local mean gap: `w = c·gap(t) = 2πc/log(t/2π)`.

## Where it provably cannot help

At `κ = gap·log X < π` this is an in-band functional —
[`band-saturation.md`](band-saturation.md)'s E4c closed that class exactly
("re-weighting primes is a choice of in-band functional, a freedom E4 already
exhausted"). So it buys nothing for `--ghy` at extreme heights, and nothing for
`--weil` (whose free functional design already spans box functionals). It is
also *not* the snap-margin experiment again: that tapered the **backward**
zero-sum (helps only off-regime, the detector being sidelobe-safe in-regime);
this tapers the **forward** prime sum, where the Gibbs tails *are* the error.

## Gate-check (Wolfram, zeros k = 1000–1100, identical brackets, truth = `ZetaZero`)

Median `|t* − γ|/gap`, plain B vs best box `c`:

| κ | X | plain | best box | gain | best c |
|---|---|---|---|---|---|
| 0.99π | 15 | 0.056 | 0.050 | 1.1× (null) | 0.375 |
| 1.49π | 59 | 0.041 | 0.016 | **2.5×** | 0.5 |
| 1.99π | 229 | 0.029 | 0.0084 | **3.5×** | 0.5 (RMS-robust 0.375) |
| 2.93π | 3000 | 0.0138 | 0.0055 | **2.5×** | 0.25 |

Both hypotheses confirmed: gains collapse to null at the band-saturation
threshold (independent consistency check of the theorem), large in the loop's
π–3π operating band; `c = 0.75` is worse everywhere (close-pair bias); optimal
`c` drifts down as κ grows. `c = 0.375` is the RMS-robust default.

## In-loop A/B (2026-07-10, 120 s wall-time each, fresh seeds, scored vs Odlyzko `zeros1`)

| run | kmin | X reached | zeros | median err/gap | p90 | primes |
|---|---|---|---|---|---|---|
| plain, anneal | 5.0→4.5 (stalled at X=104) | 8616 | 63304 | 0.0288 | 0.0688 | exact |
| box 0.375, anneal | 5.0 (never stalled) | 6886 | 57429 | **0.0127** | 0.0347 | exact |
| plain, fixed 4.0 | 4.0 | 6673 | 59028 | 0.0284 | 0.0680 | exact |
| box 0.375, fixed 4.0 | 4.0 | **7459** | 52845 | **0.0126** | 0.0346 | exact |
| box 0.375, fixed 4.5 | 4.5 | 6673 | 53242 | 0.0131 | 0.0365 | exact |

(The two anneal rows ran as one parallel pair, the three fixed rows as another;
compare within groups.) Readings:

- **2.25× forward accuracy in situ** (self-discovered primes, marching
  brackets), stable across kmin 4.0–5.0; all prime lists exact vs a sieve.
- **The accuracy converts to margin and reach.** Plain cannot hold kmin = 5.0
  (stalls immediately, anneals — over-supplying zeros to compensate noise); box
  holds it. At matched kmin = 4.0 box reaches 12% further on 10% fewer zeros
  despite the extra `sin` per term (~2× forward trig).
- Box's p90 (0.035) keeps ~2× margin below the backward snap tolerance (~0.07)
  where plain sits at the edge (0.068) — the mechanism behind the stall
  difference.

## Usage and defaults

`--loop-boxc C` sets `w = C·gap(t)`; **default 0.375 (on)**, `0` restores plain
B. The knob is runtime-only (not checkpointed): resuming an old plain-B state
continues with box-B unless `--loop-boxc 0` is passed. The ceiling message now
reports which forward step hit it.

## Follow-ups

- **Port to loopvk's forward kernel** — the weight is a per-term multiply, GPU
  trivial (unlike a Weil forward step); expected to shrink the frontier lag
  (~2–3× zeros) that dominates loopvk's cost.
- **Lower kmin floor.** With 2× extra snap margin the anneal floor (3.5 ≈ 1.1π)
  may be approachable; gains shrink toward κ = π, so expect diminishing returns,
  not a wall crossing.
- The polylog↔zeta "circularity" behind the closed form is the Möbius–ζ ladder
  (`P(s) = Σ μ(k)/k · log ζ(ks)`) — orbit §4.2.5's verdict applies: it closes
  enumeration, never evaluation; the per-prime Lerch functional equation
  degenerates to Bernoulli/Hurwitz staircases on the local (wrong) lattice.
  Theory path parked with mechanisms, as with the other in-band channels.
