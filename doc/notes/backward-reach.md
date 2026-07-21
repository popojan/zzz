# The backward step is detector-limited, not information-limited

*2026-06-21. Follow-up to [`zero-rigidity.md`](zero-rigidity.md): three negatives
closed the **forward** step (primes→zeros is information-limited below
$\kappa=\pi$ — band-saturation). But the loop's *growth rate* is set by the
**backward** step (zeros→primes), and that is a different animal. This note shows
the backward reach is **detector-limited with large unexploited slack**, so the
loop's marginality is recoverable, not a wall. Script:
`doc/ghy/reach-exponent.wls`.*

## The question

The loop map is $X_{i+1}\sim X_i^{\alpha\cdot 2\pi/\kappa}$, where $\alpha$ is the
**backward reach exponent**: from $n$ zeros (height $T$), the detector recovers
primes up to $x_{\max}\sim n^{\alpha}$. At $\kappa=\pi$ the loop is
**super-critical iff $\alpha>\tfrac12$** (the prime frontier grows without bound).
The closed-loop analysis measured $\alpha\approx0.5$–$0.6$ on the loop's own
noisy self-zeros — knife-edge. Is that a real wall or detector slack?

## Information ceiling vs measured reach

Unlike the forward step, the information here is **present**: the truncated
explicit formula
$\psi(x)=x-\sum_{|\gamma|<T}x^\rho/\rho+\dots$ has truncation error
$\sim x(\log x)^2/T$, so $n$ zeros determine prime powers up to
$x\sim T/(\log T)^2$ — **nearly linear in $n$, i.e. the ceiling is
$\alpha\to1$**. Anything below that is the detector, not missing information.

Measured (exact Odlyzko zeros, the project's envelope detector $\tau=0.30$,
flawless-prefix reach = largest $X$ with no missed prime / no composite
false-positive):

| $n$ | $\gamma_n$ | $x_{\max}$ | $x_{\max}/\gamma_n$ | local $\alpha$ |
|---|---|---|---|---|
| 500 | 811 | 309 | 0.38 | — |
| 1000 | 1419 | 599 | 0.42 | 0.96 |
| 2000 | 2515 | 1009 | 0.40 | 0.75 |
| 4000 | 4506 | 1665 | 0.37 | 0.72 |
| 8000 | 8148 | 3119 | 0.38 | 0.91 |
| 16000 | 14853 | 5231 | 0.35 | 0.75 |
| 32000 | 27260 | ≥6000 (cap) | ≥0.22 | — |

Global fit $\alpha\approx0.74$ (dragged down by the capped last point; local
slopes are $0.72$–$0.96$). So with clean zeros the **existing** detector is
already $\alpha\approx0.74$–$0.9$ — **super-critical** — and reaches
$x_{\max}\approx0.4\,\gamma_n$, far below the $\alpha=1$ ceiling. Two gaps remain:

- **$0.4 \to 1$ of $\gamma_n$**: detector slack even with perfect zeros (the
  envelope-threshold reads the raw peak; it does not deconvolve the
  Dirichlet/Fejér kernel or cancel known-prime interference).
- **$\alpha\approx0.74$ (exact) vs $0.5$–$0.6$ (loop's self-zeros)**: the cost of
  gap-scale **approximation noise** + the conservative stop. The loop's own run
  bears this out — 4450 primes ($x_{\max}=42557$) from 634k zeros
  ($\gamma\le396531$) is $x_{\max}/\gamma_n\approx0.11$, the memory's "12% of
  Nyquist," vs 0.40 here on clean zeros.

## Verdict — the real lost opportunity

The backward wall is **not information-theoretic.** It is detector- and
approximation-noise-limited, with the information ceiling at $\alpha=1$ and the
existing detector already at $\alpha\approx0.74$ on clean zeros (super-critical).
The loop's marginality ($\alpha\approx0.5$) is the **tax of running a
peak-threshold detector on its own noisy zeros**, not a barrier — and the tax is
recoverable. This is the one place in the project where the loop can beat its own
marginality (still within the fixed $\sqrt T$ *cost* per zero — efficiency vs.
direct computation is untouched; this is purely about the loop's *self-paving
rate*).

So, to "does the wall shift as you buy more approximate zeros?": the **cost** wall
(√T per zero) does not shift. But the **prime frontier** does shift, and its
growth exponent is recoverable from $\sim0.5$ toward $\sim0.74$ (clean-zero
detector) and in principle toward $1$ (information ceiling) — by a better
**backward detector** robust to gap-scale zero noise.

## Measured on the loop's OWN noisy zeros — the lever tested, and refuted

`doc/ghy/reach-detectors.wls`: flawless-prefix reach on the loop's actual
self-computed zeros (`/tmp/loop60k.txt`, the real method-B error field), three
windows, $\alpha$ = global log–log slope:

| $n$ | Dirichlet (naive) | Fejér | Jackson |
|---|---|---|---|
| 750 | 149 | 281 | 225 |
| 1500 | 464 | 461 | 419 |
| 3000 | 461 | 827 | 617 |
| 6000 | 1694 | 1543 | 827 |
| **$\alpha$** | **≈1.05** | 0.82 | 0.62 |

Two findings, both against the hypothesis above:

1. **$\alpha\approx1$ with the *naive* detector** — already super-critical on the
   loop's own noisy zeros, not the $\sim0.5$ the closed-loop analysis implied and
   not needing any clever detector. The backward step is **detector-saturated**:
   the simple full-bandwidth peak-threshold already extracts primes at a
   super-critical rate. There is **no backward lost opportunity.**
2. **Tapering *lowers* the exponent** ($1.05\to0.82\to0.62$). Mechanism: at equal
   peak height, a taper downweights the **high-$\gamma$ zeros that carry the
   Nyquist bandwidth** — exactly what reach is made of. Tapering buys *stability*
   (the Fejér/Jackson reaches are smoother, less brittle) but trades away
   reach-exponent. The snap-margin edge benefit (tolerance) does **not** transfer
   to reach.

(Caveat: 4 points, the flawless-prefix metric is brittle — a single early
sidelobe false-positive caps a reach, hence the non-monotone Dirichlet column —
so $\alpha$ is good to maybe $\pm0.2$; but $\alpha>0.5$ for the naive detector,
and the taper-hurts-exponent ordering, are robust.)

## Revised verdict

The loop's marginality is **not** backward-reach-limited. The backward direction
is super-critical and detector-saturated; the binding constraint is **forward
accuracy at the $\kappa$-frontier** — you can only make zeros accurate enough to
feed the (super-critical) backward step while $\kappa\gtrsim\pi$, and that is the
band-saturation wall ([`band-saturation.md`](band-saturation.md)), which is
information-limited and closed. So: the prime frontier *does* shift (backward
extraction is super-critical), but how far is gated by the forward $\sqrt T$ wall,
which is real and does not move except by paying $\sqrt T$ per zero. The detector
lever floated above is closed by measurement.

## Caveats

$\alpha$ measured with **exact** zeros and the flawless-prefix metric; the
$n=32000$ point hit the scan cap (so the fit underestimates); the loop's errors
are correlated common-mode, which differs from the clean case. The decisive
follow-up is $\alpha$ on the loop's *own* zeros with each candidate detector.

## Artifact index

`doc/ghy/reach-exponent.wls` (+ `.txt`). Companions:
[`zero-rigidity.md`](zero-rigidity.md), [`band-saturation.md`](band-saturation.md)
(forward = information-limited), [`snap-margin-detector.md`](snap-margin-detector.md)
(backward *tolerance*; this note is backward *reach*),
[`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md) (the loop and its
contraction-map criticality).
